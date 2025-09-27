
# -*- coding: utf-8 -*-
# Watershed Dam-System Dynamic Connectivity Builder
# -------------------------------------------------
# Requirements:
#   - geopandas, pandas, shapely, networkx
#
# Input:
#   - A Shapefile/GeoPackage of dam features (points OR polygons) with a field "淤地坝编码" (or your custom field)
#     using the hierarchical coding scheme, e.g.: 143, 14301, 1430101, 14302, ...
#     Rule: Sub-basin code = first N digits (default 3). Each hierarchical level appends two digits.
#
#   - (Optional) time fields to enable snapshots:
#       * START_FIELD: when a dam starts to "exist" (e.g., 建成年/启用年)
#       * END_FIELD  : when it stops being active (e.g., 退役年/失效年)
#     If missing, all dams are assumed to always exist.
#
# Output:
#   - A GeoPackage with layers:
#       * nodes_all, edges_all
#       * nodes_t_{t}, edges_t_{t} for each snapshot time t
#   - CSV edge lists
#
# Usage:
#   1) Edit the CONFIG block below
#   2) python watershed_dam_network.py

from typing import List, Optional, Tuple, Dict
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString

# ------------------------------
# CONFIG — EDIT THESE DEFAULTS
# ------------------------------
INPUT_PATH = r"./dams.shp"          # <-- your processed SHP (or .gpkg layer via 'path|layer')
CODE_FIELD = "淤地坝编码"               # field containing the hierarchical dam code
SUBBASIN_PREFIX_LEN = 3             # number of leading digits for the sub-basin, e.g. 143
# Optional time fields (set to None if not available)
START_FIELD = None                  # e.g. "建成年" or "start_year"
END_FIELD   = None                  # e.g. "退役年" or "end_year"
# Snapshots: list values to export dynamic states (e.g. years). Leave [] to skip.
SNAPSHOTS: List = []                # e.g. [2000, 2005, 2010, 2015, 2020]

# Output
OUT_GPKG = r"./dam_network.gpkg"
OUT_CSV_ALL = r"./dam_edges_all.csv"


# -----------------------------------------
# Helper functions for code-based topology
# -----------------------------------------
def normalize_code(v) -> str:
    """Return the code as a string without spaces; keep leading zeros if present in data."""
    if pd.isna(v):
        return ""
    return str(v).strip()


def code_level(code: str, sub_prefix_len: int) -> int:
    """Compute hierarchical level depth from code length and sub-basin prefix length.
    Each level adds 2 digits beyond the prefix.
    Example: len=3 -> level=0 (subbasin/outlet), len=5 -> level=1, len=7 -> level=2
    """
    if len(code) < sub_prefix_len:
        return -1
    tail = len(code) - sub_prefix_len
    if tail < 0 or tail % 2 != 0:
        return -1
    return tail // 2


def parent_code(code: str, sub_prefix_len: int) -> Optional[str]:
    """Parent = code with last 2 digits removed; None if already at level 0."""
    lvl = code_level(code, sub_prefix_len)
    if lvl <= 0:
        return None
    return code[:-2]


def subbasin_of(code: str, sub_prefix_len: int) -> str:
    return code[:sub_prefix_len]


# -------------------------
# Geometry / I/O utilities
# -------------------------
def read_any(path_or_layer: str) -> gpd.GeoDataFrame:
    """Read .shp or .gpkg. For GPKG specify 'path.gpkg|layername'."""
    if "|" in path_or_layer:
        path, layer = path_or_layer.split("|", 1)
        return gpd.read_file(path, layer=layer)
    else:
        return gpd.read_file(path_or_layer)


def ensure_point_geometry(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """If geometry is polygon/line, replace with centroid for connectivity drawing."""
    if gdf.geometry.geom_type.isin(["Point", "MultiPoint"]).all():
        return gdf.copy()
    out = gdf.copy()
    out["geometry"] = out.geometry.centroid
    return out


# -------------------------
# Core building functions
# -------------------------
def build_nodes_edges(
    gdf: gpd.GeoDataFrame,
    code_field: str,
    sub_prefix_len: int,
    start_field: Optional[str] = None,
    end_field: Optional[str] = None
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, pd.DataFrame]:
    """Return (nodes_all_gdf, edges_all_gdf, edges_all_df) without time filtering."""
    g = ensure_point_geometry(gdf)

    # normalize codes
    g["code"] = g[code_field].apply(normalize_code)
    g = g[g["code"] != ""].copy()

    # annotate derived attributes
    g["level"] = g["code"].apply(lambda c: code_level(c, sub_prefix_len))
    g["subbasin"] = g["code"].apply(lambda c: subbasin_of(c, sub_prefix_len))
    g["parent"] = g["code"].apply(lambda c: parent_code(c, sub_prefix_len))

    # build a quick lookup for parent geometry
    geom_map: Dict[str, Point] = dict(zip(g["code"], g.geometry))

    # construct edges as lines from child to parent (if parent exists in dataset)
    edge_records = []
    for row in g.itertuples(index=False):
        if row.parent is None:
            continue
        if row.parent in geom_map:
            p1 = row.geometry
            p0 = geom_map[row.parent]
            line = LineString([p1, p0])
            edge_records.append({
                "child": row.code,
                "parent": row.parent,
                "level": row.level,
                "subbasin": row.subbasin,
                "geometry": line
            })
        else:
            # parent not present in data; still record attribute-only edge
            edge_records.append({
                "child": row.code,
                "parent": row.parent,
                "level": row.level,
                "subbasin": row.subbasin,
                "geometry": None
            })

    edges_gdf = gpd.GeoDataFrame(edge_records, geometry="geometry", crs=g.crs)
    edges_df = edges_gdf.drop(columns=["geometry"]).copy()

    # Keep time columns if provided
    keep_cols = ["code", "level", "subbasin", "parent"]
    if start_field and start_field in g.columns:
        keep_cols.append(start_field)
    if end_field and end_field in g.columns:
        keep_cols.append(end_field)
    nodes = g[keep_cols + ["geometry"]].copy()

    return nodes, edges_gdf, edges_df


def filter_by_time(
    nodes: gpd.GeoDataFrame,
    edges: gpd.GeoDataFrame,
    t
) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Filter nodes/edges to those 'active' at time t using START_FIELD/END_FIELD if present."""
    n = nodes.copy()
    if START_FIELD and START_FIELD in n.columns:
        n = n[(n[START_FIELD].isna()) | (n[START_FIELD] <= t)]
    if END_FIELD and END_FIELD in n.columns:
        n = n[(n[END_FIELD].isna()) | (n[END_FIELD] > t)]
    active_codes = set(n["code"])
    e = edges[edges["child"].isin(active_codes) & edges["parent"].isin(active_codes)].copy()
    return n, e


def export_all(nodes_all, edges_all, out_gpkg: str, out_csv_all: str, snapshots: List):
    """Save everything. For snapshots, create per-t layers/CSVs."""
    out_path = Path(out_gpkg)
    if out_path.exists():
        out_path.unlink()  # overwrite cleanly

    nodes_all.to_file(out_gpkg, layer="nodes_all", driver="GPKG")
    edges_all.to_file(out_gpkg, layer="edges_all", driver="GPKG")
    edges_all.drop(columns=["geometry"]).to_csv(out_csv_all, index=False, encoding="utf-8-sig")

    # snapshots
    for t in snapshots:
        n_t, e_t = filter_by_time(nodes_all, edges_all, t)
        n_t.to_file(out_gpkg, layer=f"nodes_t_{t}", driver="GPKG")
        e_t.to_file(out_gpkg, layer=f"edges_t_{t}", driver="GPKG")
        e_t.drop(columns=["geometry"]).to_csv(out_path.with_name(f"dam_edges_t_{t}.csv"), index=False, encoding="utf-8-sig")


def main():
    gdf = read_any(INPUT_PATH)
    nodes_all, edges_all, _ = build_nodes_edges(
        gdf, CODE_FIELD, SUBBASIN_PREFIX_LEN, START_FIELD, END_FIELD
    )
    export_all(nodes_all, edges_all, OUT_GPKG, OUT_CSV_ALL, SNAPSHOTS)
    print(f"Done. Wrote: {OUT_GPKG} and {OUT_CSV_ALL}")
    if SNAPSHOTS:
        print(f"Snapshots: {SNAPSHOTS}")


if __name__ == "__main__":
    main()
