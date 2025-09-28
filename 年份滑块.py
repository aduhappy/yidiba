# -*- coding: utf-8 -*-
# 年份滑块 + 子流域（弹窗仅显示：坝编码、编码、淤地坝）
import re
import argparse
import geopandas as gpd
gpd.options.io_engine = "fiona"
import folium
from folium.plugins import MarkerCluster
from folium.features import DivIcon
from branca.element import Element

# ===== 默认路径（按需改成你的实际路径）=====
BASE = r"F:\OneDrive - 8f7hnp\work\2025年\9月\孙彭成老师\To Lintao\python\淤地坝"
DEFAULT_GPKG   = BASE + r"\dam_network.gpkg"
DEFAULT_WS_SHP = BASE + r"\Watershed.shp"
DEFAULT_ATTRS  = BASE + r"\dams_clean.gpkg"
DEFAULT_HTML   = BASE + r"\dam_network_timeslider_two_selects_nameonly.html"
# =======================================

def normalize_columns(df):
    m = {}
    for c in df.columns:
        cc = str(c)
        cc2 = (cc.replace("\ufeff","").replace("\xa0","")
                 .replace(" ","").replace("\u3000","").strip())
        m[c] = cc2
    if any(k!=v for k,v in m.items()):
        df = df.rename(columns=m, inplace=False)
    return df

def resolve_column_name(df, name):
    if name is None: return None
    cols = list(df.columns)
    if name in cols: return name
    nm = (str(name).replace("\ufeff","").replace("\xa0","")
                    .replace(" ","").replace("\u3000","").strip())
    norm = { (str(c).replace("\ufeff","").replace("\xa0","")
                   .replace(" ","").replace("\u3000","").strip()): c for c in cols }
    return norm.get(nm, None)

def coerce_year(v):
    if v is None: return None
    s = str(v)
    m = re.search(r"(18|19|20)\d{2}", s)
    if not m: return None
    y = int(m.group(0))
    return y if 1800 <= y <= 2100 else None

def pick_year_column(df, prefer_list=("建成年","建成年份","年份","年","竣工年","修建年")):
    cols = list(df.columns)
    for pref in prefer_list:
        c = resolve_column_name(df, pref)
        if c: return c
    cands = [c for c in cols if any(p in str(c) for p in prefer_list)]
    def score(c):
        try: return int(df[c].map(coerce_year).notna().sum())
        except Exception: return 0
    if cands:
        s = sorted([(score(c), c) for c in cands], reverse=True)
        if s and s[0][0] > 0: return s[0][1]
    scores = []
    for c in cols:
        try: scores.append((int(df[c].map(coerce_year).notna().sum()), c))
        except Exception: pass
    scores.sort(reverse=True)
    if scores and scores[0][0] > 0: return scores[0][1]
    return None

def find_code_col(df):
    for c in ["code","编码","淤地坝编码","坝码","id","ID","Code","CODE"]:
        if c in df.columns: return c
    return None

def prepare_code_series(s):
    return (s.astype(str)
              .str.replace("\ufeff","", regex=False)
              .str.replace("\xa0","", regex=False)
              .str.replace("\u3000","", regex=False)
              .str.strip()
              .str.replace(".0","", regex=False))

def attach_year(nodes, prefer_col=None, attrs_gpkg=None, attrs_layer="dams_points"):
    """
    仅补：__year__（建成年）、ext_code（外部长码“编码”）、dam_name（“淤地坝”名称）
    """
    def _norm(df): return normalize_columns(df)
    def _pick(df, cands):
        for k in cands:
            if k in df.columns: return k
        for k in df.columns:
            sk = str(k).replace("\ufeff","").replace("\xa0","").replace("\u3000","").replace(" ","")
            for want in cands:
                if want.replace(" ","") in sk: return k
        return None

    nodes = _norm(nodes)

    # nodes 自身建成年
    prefer_col = resolve_column_name(nodes, prefer_col)
    ycol_nodes = prefer_col or pick_year_column(nodes)
    if ycol_nodes:
        yrs = nodes[ycol_nodes].map(coerce_year)
        if yrs.notna().sum() > 0:
            nodes["__year__"] = yrs
            print(f"[attach_year] use nodes.'{ycol_nodes}'  non-null={int(yrs.notna().sum())}/{len(nodes)}")
        else:
            print(f"[attach_year] nodes.'{ycol_nodes}' exists but all invalid/empty")

    # 属性库合并（按“坝编码”对齐）
    if attrs_gpkg:
        try:
            attrs = gpd.read_file(attrs_gpkg, layer=attrs_layer)
        except Exception:
            try: attrs = gpd.read_file(attrs_gpkg)
            except Exception: attrs = None
        if attrs is not None:
            attrs = _norm(attrs)
            cn = find_code_col(nodes)
            ca = find_code_col(attrs)
            ay = pick_year_column(attrs)  # 年份
            a_ext  = _pick(attrs, ["编码","统一编码","全码","唯一编码"])
            a_name = _pick(attrs, ["淤地坝","名称","坝名","工程名"])  # 仅此做弹窗

            if cn and ca:
                keep = [c for c in [ca, ay, a_ext, a_name] if c]
                tmp = attrs[keep].copy()
                left = prepare_code_series(nodes[cn])
                tmp[ca] = prepare_code_series(tmp[ca])
                if ay: tmp[ay] = tmp[ay].map(coerce_year)
                before = len(nodes)
                nodes = nodes.merge(tmp, how="left", left_on=left, right_on=tmp[ca])
                if ay:
        if ay in nodes.columns:
            nodes["__year__"] = nodes["__year__"].fillna(nodes[ay])
        else:
            candidates = [c for c in nodes.columns if str(c).startswith(str(ay))]
            if candidates:
                print(f"[warn] year_col '{ay}' not found, use '{candidates[0]}' instead")
                nodes["__year__"] = nodes["__year__"].fillna(nodes[candidates[0]])
                if a_ext:  nodes["ext_code"] = nodes[a_ext]
                if a_name: nodes["dam_name"] = nodes[a_name]
                ok = int(nodes["__year__"].notna().sum()) if "__year__" in nodes.columns else 0
                print(f"[attach_year] merge attrs by nodes.{cn} ~ attrs.{ca} | year='{ay}' non-null={ok}/{before}")
            else:
                print("[attach_year] cannot find code field in nodes/attrs")
        else:
            print("[attach_year] read attrs failed:", attrs_gpkg)

    if "__year__"  not in nodes.columns: nodes["__year__"]  = None
    if "ext_code"  not in nodes.columns: nodes["ext_code"]  = None
    if "dam_name"  not in nodes.columns: nodes["dam_name"]  = None
    return nodes

def pick_ws_field(df):
    df = normalize_columns(df)
    for cand in ["name","Name","NAME","ws_name","WS_NAME","Subbasin","SUBBASIN","HUC","HUC12","code","CODE","编号","子流域","小流域"]:
        if cand in df.columns: return cand
    for c in df.columns:
        if df[c].dtype == "object": return c
    return None

def level_color(val):
    try:
        if isinstance(val, dict):
            for k in ("value", "level"):
                if k in val:
                    val = val[k]; break
        val = int(val)
    except Exception:
        val = None
    return {0:"green",1:"orange",2:"lightblue",3:"purple"}.get(val, "red")

def add_all_points_debug(m, nodes_ws):
    fg = folium.FeatureGroup(name="所有坝点(调试)", show=True, control=True)
    mc = MarkerCluster(name="cluster").add_to(fg)
    for r in nodes_ws.itertuples(index=False):
        p = getattr(r, "geometry", None)
        if p is None or p.is_empty: continue
        pt = p if p.geom_type == "Point" else p.centroid
        lat, lon = pt.y, pt.x
        code = getattr(r,"code","")
        popup = (
            f"<b>坝编码</b>: {code}<br>"
            f"<b>编码</b>: {getattr(r,'ext_code','') or ''}<br>"
            f"<b>淤地坝</b>: {getattr(r,'dam_name','') or ''}"
        )
        folium.CircleMarker(
            (lat, lon), radius=4, color="#ff7800", fill=True, fill_opacity=0.7,
            popup=popup, tooltip=folium.Tooltip(f"坝编码:{code} | 编码:{getattr(r,'ext_code','') or '-'}")
        ).add_to(mc)
    fg.add_to(m)
    return fg

def add_code_labels_layer(m, gdf, name="编码标注"):
    fg = folium.FeatureGroup(name=name, show=False, control=True)
    for r in gdf.itertuples(index=False):
        p = getattr(r, "geometry", None)
        if p is None or p.is_empty: continue
        pt = p if p.geom_type == "Point" else p.centroid
        code = getattr(r,"code","") or ""
        ext  = getattr(r,"ext_code","") or ""
        if not code and not ext: continue
        html = f'<div style="font-size:10px;color:#222;background:rgba(255,255,255,.75);padding:1px 2px;border:1px solid #888;border-radius:2px;">{code}|{ext}</div>'
        folium.Marker((pt.y, pt.x), icon=DivIcon(html=html)).add_to(fg)
    fg.add_to(m); return fg

def build_group(m, nodes_ws, edges, year, ws_name):
    fg = folium.FeatureGroup(name=f"≤{year} | {ws_name}", show=False, control=False)
    nn = nodes_ws.copy()
    if nn["__year__"].notna().sum() > 0:
        nn = nn[nn["__year__"] <= year]
    if ws_name != "全部":
        nn = nn[nn["__ws__"] == ws_name]
    if "child" in edges.columns and "parent" in edges.columns:
        active = set(nn["code"]) if "code" in nn.columns else set(nn.index.astype(str))
        ee = edges[edges["child"].isin(active) & edges["parent"].isin(active)]
        for r in ee.itertuples(index=False):
            g = getattr(r, "geometry", None)
            if g is None or g.is_empty: continue
            if g.geom_type == "MultiLineString":
                for ls in g.geoms:
                    folium.PolyLine([(y,x) for x,y in ls.coords], weight=2, opacity=0.8).add_to(fg)
            elif g.geom_type == "LineString":
                folium.PolyLine([(y,x) for x,y in g.coords], weight=2, opacity=0.8).add_to(fg)
    for r in nn.itertuples(index=False):
        p = getattr(r, "geometry", None)
        if p is None or p.is_empty: continue
        pt = p if p.geom_type == "Point" else p.centroid
        lat, lon = pt.y, pt.x
        code = getattr(r,"code","")
        popup = (
            f"<b>坝编码</b>: {code}<br>"
            f"<b>编码</b>: {getattr(r,'ext_code','') or ''}<br>"
            f"<b>淤地坝</b>: {getattr(r,'dam_name','') or ''}"
        )
        folium.CircleMarker(
            (lat, lon), radius=5, color=level_color(getattr(r,"level",None)),
            fill=True, fill_opacity=0.9, popup=popup,
            tooltip=folium.Tooltip(f"坝编码:{code} | 编码:{getattr(r,'ext_code','') or '-'}")
        ).add_to(fg)
    fg.add_to(m); return fg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gpkg",       default=DEFAULT_GPKG)
    ap.add_argument("--ws-shp",     default=DEFAULT_WS_SHP)
    ap.add_argument("--out",        default=DEFAULT_HTML)
    ap.add_argument("--year-col",   default=None)
    ap.add_argument("--years",      nargs="*")
    ap.add_argument("--ws-field",   default=None)
    ap.add_argument("--ws-max",     type=int, default=30)
    ap.add_argument("--attrs-gpkg", default=DEFAULT_ATTRS)
    ap.add_argument("--attrs-layer",default="dams_points")
    args = ap.parse_args()

    nodes = gpd.read_file(args.gpkg, layer="nodes_all")
    edges = gpd.read_file(args.gpkg, layer="edges_all")
    ws    = gpd.read_file(args.ws_shp)

    for g in (nodes, edges, ws):
        try:
            if g.crs is not None and g.crs.to_epsg() != 4326: g.to_crs(4326, inplace=True)
        except Exception:
            g.to_crs(4326, inplace=True)

    nodes = attach_year(nodes, prefer_col=args.year_col, attrs_gpkg=args.attrs_gpkg, attrs_layer=args.attrs_layer)

    wcol = args.ws_field or pick_ws_field(ws)
    ws["__ws__"] = (ws[wcol].astype(str) if wcol else ws.index.astype(str))

    try:
        nodes_ws = gpd.sjoin(nodes, ws[["__ws__","geometry"]], how="left", predicate="intersects")
    except Exception:
        nodes_ws = gpd.sjoin(nodes, ws[["__ws__","geometry"]], how="left")
    nodes_ws["__ws__"] = nodes_ws["__ws__"].fillna("未分配")

    if args.years:
        years = []
        for token in args.years:
            for p in str(token).replace("，", ",").split(","):
                p = p.strip()
                if p: years.append(int(p))
        years = sorted(set(years))
    else:
        years = sorted(set([int(y) for y in nodes_ws["__year__"].dropna().tolist()]))
        if not years: years = [2000]
        elif len(years) > 25:
            ymin, ymax = min(years), max(years)
            step = max(1, (ymax - ymin)//18)
            years = list(range(ymin, ymax+1, step))

    top = nodes_ws.groupby("__ws__")["geometry"].count().sort_values(ascending=False).head(max(1, int(args.ws_max))).index.tolist()
    ws_list = ["全部"] + top

    try:
        c = (nodes_ws.geometry.union_all() if hasattr(nodes_ws.geometry,"union_all") else nodes_ws.geometry.unary_union).centroid
        center = [c.y, c.x]
    except Exception:
        center = [35.0,109.0]
    m = folium.Map(location=center, zoom_start=9, tiles=None)
    folium.TileLayer(tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
                     attr="© Esri", name="Esri World Imagery").add_to(m)

    ws_layer = folium.GeoJson(ws, name="子流域边界",
                              style_function=lambda feat: {"color":"#888","weight":1,"fill":False})
    ws_layer.add_to(m); ws_var = ws_layer.get_name()

    add_all_points_debug(m, nodes_ws)
    add_code_labels_layer(m, nodes_ws, name="编码标注")

    layer_map = {}
    for y in years:
        layer_map[str(y)] = {}
        for w in ws_list:
            nm = build_group(m, nodes_ws, edges, y, w).get_name()
            layer_map[str(y)][w] = nm

    def js_obj(d):
        lines = []
        for y, sub in d.items():
            inner = ", ".join(['"%s": %s' % (w, sub[w]) for w in sub])
            lines.append('  "%s": { %s }' % (y, inner))
        return "{\n" + ",\n".join(lines) + "\n}"
    layer_js = js_obj(layer_map)
    year_list_js = "[" + ",".join([str(y) for y in years]) + "]"

    ws_opts = "".join(['<option value="%s">%s</option>' % (w, w) for w in ws_list])
    panel_html = f'''<div id="wsPanel" style="position:absolute;right:12px;top:12px;z-index:1000;background:#fff;padding:8px 10px;
box-shadow:0 1px 5px rgba(0,0,0,.4);border-radius:6px;font:12px/1.4 sans-serif;">
  <div style="margin-bottom:4px;">子流域</div>
  <select id="wsSel" style="padding:2px 4px;width:160px;margin-bottom:6px;">{ws_opts}</select>
  <div style="margin-top:6px;display:flex;gap:10px;align-items:center;">
    <label style="display:flex;gap:4px;align-items:center;"><input id="btnWS" type="checkbox" checked> 子流域边界</label>
  </div>
</div>
<div id="yrPanel" style="position:absolute;left:12px;bottom:12px;z-index:1000;background:#fff;padding:6px 10px;
box-shadow:0 1px 5px rgba(0,0,0,.4);border-radius:6px;font:12px/1.4 sans-serif;">
  <div id="yearLabel" style="text-align:center;margin-bottom:6px;"></div>
  <input id="yearSlider" type="range" min="0" max="{len(years)-1}" step="1" value="{len(years)-1}" style="width:260px;">
</div>'''
    m.get_root().html.add_child(Element(panel_html))

    js_code = '''
(function(){
  var LAYERS = __LAYER_JS__;
  var YEARS  = __YEAR_LIST__;
  function getMap(){
    for (var k in window){ try{ if (window[k] && window[k] instanceof L.Map) return window[k]; }catch(e){} }
    return null;
  }
  var map = getMap(); if (!map) return;
  var wsSel   = document.getElementById('wsSel');
  var slider  = document.getElementById('yearSlider');
  var btnWS   = document.getElementById('btnWS');
  var yearLab = document.getElementById('yearLabel');
  var active  = null;
  function addL(x){ try{ map.addLayer(x); }catch(e){} }
  function rmL(x){ try{ map.removeLayer(x); }catch(e){} }
  function update(){
    var idx = parseInt(slider.value,10), y = YEARS[idx], ws = wsSel.value;
    yearLab.innerHTML = '&le; ' + y + ' 年';
    var next = (LAYERS[String(y)]||{})[ws] || null;
    if (active) rmL(active);
    active = next;
    if (active) addL(active);
  }
  btnWS.addEventListener('change', function(){ if (this.checked) addL(__WS_VAR__); else rmL(__WS_VAR__); });
  wsSel.addEventListener('change', update);
  slider.addEventListener('input', update);
  addL(__WS_VAR__);
  update();
})();
'''
    js_code = js_code.replace("__LAYER_JS__", layer_js).replace("__YEAR_LIST__", year_list_js).replace("__WS_VAR__", ws_var)
    m.get_root().script.add_child(Element(js_code))

    folium.LayerControl(collapsed=False).add_to(m)
    m.save(args.out)
    print("[done] Saved:", args.out)

if __name__ == "__main__":
    main()
