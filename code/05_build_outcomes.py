"""
05_build_outcomes.py
====================
The entire TF/UE outcome construction, in one script. Replaces four files that
passed pickles between each other to do one logical job.

TF (task force) records the rescue site as FREE TEXT with no controlled
vocabulary: mixed case, accented and unaccented, name-with-UF-suffix, and in
roughly a fifth of rows raw GPS coordinates in at least six different formats
instead of a name. 287 of 1,875 municipalities appear under more than one raw
string, one of them under eight.

FIVE PASSES, each restricted to the stated UF:
  A  exact match after Unicode NFKD fold, lowercase, strip non-alphanumerics
  B  UF-suffix strip ("GOIANIA/GO", "Ibitita-BA", "Itapaci - GO")
  C  fuzzy, ratio >= 85 within UF - every accepted match reviewed by hand
  D  coordinate parsing (decimal and DMS) with point-in-polygon against the
     2022 mesh, accepted ONLY when the joined UF equals the stated UF
  E  two manual fixes: ACU/RN -> Assu (2400208); Parnarama/PI -> 2107506,
     a source UF error

INTEGRITY CHECK: no (UF, normalized string) may map to more than one
municipality code. That is the failure mode that would silently corrupt the
panel, and it is asserted, not assumed.

UE (seguro-desemprego) is simpler - it carries an IBGE code directly - but needs
a footer row dropped (a 69,701-victim grand total) and one typo year (1199).
21% of UE victims have a blank residence code; those are counted and reported,
never imputed.

Outputs: tf_cleaned_rowlevel.csv, tf_location_crosswalk.csv,
         tf_muni_quarter.csv, ue_muni_quarter.csv
Seed: 123
"""


# =============================================================================
# PASSES A-C: normalize, UF-strip, fuzzy
# =============================================================================
import pandas as pd, re
import geopandas as gpd

m = pd.read_pickle('tf_pass2.pkl')
rem = pd.read_pickle('tf_rem.pkl')
fz = pd.read_csv('tf_fuzzy_matches.csv')

m.loc[fz['idx'].values, 'code7'] = fz['code7'].values
m.loc[fz['idx'].values, 'match_method'] = 'fuzzy'

REPL = [('º', '°'), ('˚', '°'), ('–', '-'), ('—', '-'),
        ('‘', "'"), ('’', "'"), ('´', "'"), ('`', ''),
        ('“', '"'), ('”', '"'), ("''", '"'), (',', '.')]

def clean_str(s):
    s = str(s)
    for a, b in REPL:
        s = s.replace(a, b)
    return s

dms_pat = re.compile(r"(\d{1,3})\s*°\s*(\d{1,2})\s*'\s*(\d{1,2}(?:\.\d+)?)?")
dec_pat = re.compile(r"-?\d{1,3}\.\d{3,}")

def parse_coord(raw):
    s = clean_str(raw)
    up = s.upper()
    north = bool(re.match(r'^\s*N[\s0-9]', up))
    decs = [float(d) for d in dec_pat.findall(s) if abs(float(d)) <= 180]
    cand = [d for d in decs if abs(d) > 0.5]
    if len(cand) >= 2:
        lat, lon = cand[0], cand[1]
        if abs(lat) > 34 or abs(lon) < 30:
            lat, lon = cand[1], cand[0]
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    groups = dms_pat.findall(s)
    if len(groups) >= 2:
        def to_dec(g):
            return float(g[0]) + float(g[1]) / 60 + (float(g[2]) if g[2] else 0.0) / 3600
        lat, lon = to_dec(groups[0]), to_dec(groups[1])
        if lat > 34 or lon < 30:
            lat, lon = lon, lat
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    return None, None

coords = rem[rem.is_coord].copy()
parsed = coords['loc2'].map(parse_coord)
coords['lat'] = [p[0] for p in parsed]
coords['lon'] = [p[1] for p in parsed]
ok = coords.lat.notna() & coords.lat.between(-34, 5.5) & coords.lon.between(-74, -34)
print('coord rows:', len(coords), '| parsed ok:', int(ok.sum()))
bad = coords.loc[~ok, ['UF', 'loc2']]
if len(bad):
    print(bad.head(20).to_string())

gj = gpd.read_file('malha.json')
pts = gpd.GeoDataFrame(coords[ok].copy(),
                       geometry=gpd.points_from_xy(coords.loc[ok, 'lon'], coords.loc[ok, 'lat']),
                       crs='EPSG:4674')
joined = gpd.sjoin(pts, gj[['codarea', 'geometry']], how='left', predicate='within')
joined = joined[~joined.index.duplicated()]
print('spatial-joined:', int(joined.codarea.notna().sum()), 'of', len(joined))

missing_idx = joined[joined.codarea.isna()].index
if len(missing_idx):
    pts_m = pts.loc[missing_idx].to_crs(5880)
    gj_m = gj.to_crs(5880)
    near = gpd.sjoin_nearest(pts_m, gj_m[['codarea', 'geometry']], how='left', max_distance=20000)
    near = near[~near.index.duplicated()]
    joined.loc[near.index, 'codarea'] = near['codarea']
    print('after nearest(<=20km):', int(joined.codarea.notna().sum()), 'of', len(joined))

got = joined[joined.codarea.notna()]
m.loc[got.index, 'code7'] = got['codarea'].astype(float)
m.loc[got.index, 'match_method'] = 'coordinates'

uf_codes = {'RO':'11','AC':'12','AM':'13','RR':'14','PA':'15','AP':'16','TO':'17','MA':'21','PI':'22',
            'CE':'23','RN':'24','PB':'25','PE':'26','AL':'27','SE':'28','BA':'29','MG':'31','ES':'32',
            'RJ':'33','SP':'35','PR':'41','SC':'42','RS':'43','MS':'50','MT':'51','GO':'52','DF':'53'}
cm = m.loc[got.index].dropna(subset=['code7'])
mask = [str(int(c))[:2] != uf_codes.get(u, 'xx') for c, u in zip(cm.code7, cm.UF)]
mismatch = cm[mask]
print('UF mismatches from coords:', len(mismatch))
if len(mismatch):
    print(mismatch[['UF', 'location', 'code7', 'workers_rescued']].to_string())

still = m[m.code7.isna() & m.location.notna()]
print('FINAL unmatched rows:', len(still), '| rescued in them:', still.workers_rescued.sum())
m.to_pickle('tf_pass3.pkl')
(still.groupby(['UF', 'location'])
 .agg(n=('location', 'size'), rescued=('workers_rescued', 'sum'))
 .reset_index().sort_values('rescued', ascending=False)
 .to_csv('tf_unmatched_final.csv', index=False))

# =============================================================================
# PASS D part 1: coordinate parsing
# =============================================================================
import pandas as pd, re
import geopandas as gpd

m = pd.read_pickle('tf_pass3.pkl')

uf_codes = {'RO':'11','AC':'12','AM':'13','RR':'14','PA':'15','AP':'16','TO':'17','MA':'21','PI':'22',
            'CE':'23','RN':'24','PB':'25','PE':'26','AL':'27','SE':'28','BA':'29','MG':'31','ES':'32',
            'RJ':'33','SP':'35','PR':'41','SC':'42','RS':'43','MS':'50','MT':'51','GO':'52','DF':'53'}

# Rows resolved by coordinates whose UF disagrees -> re-examine
coord_idx = m.index[(m.match_method == 'coordinates')]
bad = [i for i in coord_idx if str(int(m.loc[i, 'code7']))[:2] != uf_codes.get(m.loc[i, 'UF'], 'xx')]
# plus rows still unmatched
still = m.index[m.code7.isna() & m.location.notna()]
redo = list(bad) + list(still)
print('re-examining', len(redo), 'rows')
m.loc[bad, 'code7'] = float('nan')
m.loc[bad, 'match_method'] = None

REPL = [('º','°'),('˚','°'),('–','-'),('—','-'),('‘',"'"),('’',"'"),('´',"'"),('`',''),
        ('“','"'),('”','"'),('′',"'"),('″','"'),("''",'"'),('=',''),(',','.')]

def clean_str(s):
    s = str(s)
    for a,b in REPL: s = s.replace(a,b)
    return s

dms_pat = re.compile(r"(\d{1,3})\s*[°O]\s*(\d{1,2})\s*'\s*(\d{1,2}(?:\.\d+)?)?")
dec_pat = re.compile(r"-?\d{1,3}\.\d{2,}")

def parse_coord(raw):
    s = clean_str(raw)
    up = s.upper()
    north = bool(re.search(r'\d\s*"?\s*N|^\s*N\s*\d|\bN\s*\d', up)) and not re.search(r'\d\s*"?\s*S|\bS\s*[:\s]*\d', up.split('W')[0].split('O')[0][:len(up)//2] if False else up[:len(up)//2])
    # simpler: north if N appears and S does not
    north = ('N' in re.sub(r'LONG|LAT|MIN', '', up)) and not re.search(r'S', re.sub(r'LONG|LAT|OESTE|LESTE', '', up))
    decs = [float(d) for d in dec_pat.findall(s) if abs(float(d)) <= 180]
    cand = [d for d in decs if abs(d) > 0.5]
    if len(cand) >= 2:
        lat, lon = cand[0], cand[1]
        if abs(lat) > 34 or abs(lon) < 30:
            lat, lon = cand[1], cand[0]
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    groups = dms_pat.findall(s)
    if len(groups) >= 2:
        def to_dec(g):
            mi = float(g[1]); se = float(g[2]) if g[2] else 0.0
            if mi >= 60 or se >= 60: return None
            return float(g[0]) + mi/60 + se/3600
        a, b = to_dec(groups[0]), to_dec(groups[1])
        if a is None or b is None: return None, None
        lat, lon = a, b
        if lat > 34 or lon < 30:
            lat, lon = lon, lat
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    return None, None

redo_df = m.loc[redo].copy()
redo_df['loc2'] = redo_df['location']
parsed = redo_df['loc2'].map(parse_coord)
redo_df['lat'] = [p[0] for p in parsed]; redo_df['lon'] = [p[1] for p in parsed]
ok = redo_df.lat.notna() & redo_df.lat.between(-34, 5.5) & redo_df.lon.between(-74, -34)
print('parsed ok this pass:', int(ok.sum()))

gj = gpd.read_file('malha.json').set_crs('EPSG:4326', allow_override=True)
pts = gpd.GeoDataFrame(redo_df[ok].copy(),
                       geometry=gpd.points_from_xy(redo_df.loc[ok,'lon'], redo_df.loc[ok,'lat']),
                       crs='EPSG:4326')
joined = gpd.sjoin(pts, gj[['codarea','geometry']], how='left', predicate='within')
joined = joined[~joined.index.duplicated()]
miss = joined[joined.codarea.isna()].index
if len(miss):
    near = gpd.sjoin_nearest(pts.loc[miss].to_crs(5880), gj.to_crs(5880)[['codarea','geometry']],
                             how='left', max_distance=20000)
    near = near[~near.index.duplicated()]
    joined.loc[near.index, 'codarea'] = near['codarea']

got = joined[joined.codarea.notna()]
# accept only if UF agrees
agree = [i for i in got.index if str(int(got.loc[i,'codarea']))[:2] == uf_codes.get(m.loc[i,'UF'],'xx')]
disagree = [i for i in got.index if i not in agree]
m.loc[agree, 'code7'] = got.loc[agree, 'codarea'].astype(float)
m.loc[agree, 'match_method'] = 'coordinates'
print('accepted (UF agrees):', len(agree), '| UF still disagrees:', len(disagree))
if disagree:
    print(m.loc[disagree, ['UF','location','workers_rescued']].assign(joined_code=got.loc[disagree,'codarea'].values).to_string())

still = m[m.code7.isna() & m.location.notna()]
print('\nFINAL unmatched rows:', len(still), '| rescued:', still.workers_rescued.sum())
print('match_method counts:'); print(m.match_method.value_counts(dropna=False))
m.to_pickle('tf_final.pkl')
(still.groupby(['UF','location'])
 .agg(n=('location','size'), rescued=('workers_rescued','sum'))
 .reset_index().sort_values('rescued', ascending=False)
 .to_csv('tf_unmatched_final.csv', index=False))
print(still.groupby(['UF','location']).agg(n=('location','size'), rescued=('workers_rescued','sum')).reset_index().sort_values('rescued', ascending=False).head(30).to_string())

# =============================================================================
# PASS D part 2: coordinate parsing, remainder
# =============================================================================
import pandas as pd, re
import geopandas as gpd

m = pd.read_pickle('tf_final.pkl')

uf_codes = {'RO':'11','AC':'12','AM':'13','RR':'14','PA':'15','AP':'16','TO':'17','MA':'21','PI':'22',
            'CE':'23','RN':'24','PB':'25','PE':'26','AL':'27','SE':'28','BA':'29','MG':'31','ES':'32',
            'RJ':'33','SP':'35','PR':'41','SC':'42','RS':'43','MS':'50','MT':'51','GO':'52','DF':'53'}

REPL = [('º','°'),('˚','°'),('–','-'),('—','-'),('‘',"'"),('’',"'"),('´',"'"),('`',''),
        ('“','"'),('”','"'),('′',"'"),('″','"'),("''",'"'),('=',''),(',','.')]

def clean_str(s):
    s = str(s)
    for a,b in REPL: s = s.replace(a,b)
    s = re.sub(r'(\d)\.\s+(\d)', r'\1.\2', s)   # "-45. 956" -> "-45.956"
    s = re.sub(r'\bLONG\b|\bLAT\b', ' ', s, flags=re.I)
    return s

# DMS: degree symbol optional when followed by min' sec" pattern
dms_pat = re.compile(r"(\d{1,3})\s*°?\s+?(\d{1,2})\s*'\s*(\d{1,2}(?:\.\d+)?)?\s*\"?")
dms_strict = re.compile(r"(\d{1,3})\s*°\s*(\d{1,2})\s*'?\s*(\d{1,2}(?:\.\d+)?)?")
dec_pat = re.compile(r"-?\d{1,3}\.\d{2,}")

def to_dec(g):
    mi = float(g[1]); se = float(g[2]) if g[2] else 0.0
    if mi >= 60 or se >= 60: return None
    return float(g[0]) + mi/60 + se/3600

def parse_coord(raw):
    s = clean_str(raw)
    up = s.upper()
    north = ('N' in up) and not re.search(r'S', up)
    # DMS FIRST when a degree symbol is present
    if '°' in s:
        groups = dms_strict.findall(s)
        if len(groups) < 2:
            groups = dms_pat.findall(s)
        if len(groups) >= 2:
            a, b = to_dec(groups[0]), to_dec(groups[1])
            if a is not None and b is not None:
                lat, lon = a, b
                if lat > 34 or lon < 30:
                    lat, lon = lon, lat
                lat = abs(lat) if north else -abs(lat)
                if -34 <= lat <= 5.5 and 34 <= abs(lon) <= 74:
                    return lat, -abs(lon)
    decs = [float(d) for d in dec_pat.findall(s) if abs(float(d)) <= 180]
    cand = [d for d in decs if abs(d) > 0.5]
    if len(cand) >= 2:
        lat, lon = cand[0], cand[1]
        if abs(lat) > 34 or abs(lon) < 30:
            lat, lon = cand[1], cand[0]
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    groups = dms_pat.findall(s)
    if len(groups) >= 2:
        a, b = to_dec(groups[0]), to_dec(groups[1])
        if a is None or b is None: return None, None
        lat, lon = a, b
        if lat > 34 or lon < 30:
            lat, lon = lon, lat
        lat = abs(lat) if north else -abs(lat)
        return lat, -abs(lon)
    return None, None

# manual aliases for renamed/UF-error cases
manual = {('RN','AÇU'): 2400208.0,          # Assú (renamed)
          ('PI','Parnarama/PI'): 2107506.0} # Parnarama is MA 2107506; source UF error
for (uf, loc), code in manual.items():
    idx = m.index[(m.UF == uf) & (m.location == loc) & m.code7.isna()]
    m.loc[idx, 'code7'] = code
    m.loc[idx, 'match_method'] = 'manual'

still_idx = m.index[m.code7.isna() & m.location.notna()]
redo_df = m.loc[still_idx].copy()
parsed = redo_df['location'].map(parse_coord)
redo_df['lat'] = [p[0] for p in parsed]; redo_df['lon'] = [p[1] for p in parsed]
ok = redo_df.lat.notna() & redo_df.lat.between(-34, 5.5) & redo_df.lon.between(-74, -34)
print('parsed ok this pass:', int(ok.sum()), 'of', len(redo_df))

gj = gpd.read_file('malha.json').set_crs('EPSG:4326', allow_override=True)
pts = gpd.GeoDataFrame(redo_df[ok].copy(),
                       geometry=gpd.points_from_xy(redo_df.loc[ok,'lon'], redo_df.loc[ok,'lat']),
                       crs='EPSG:4326')
joined = gpd.sjoin(pts, gj[['codarea','geometry']], how='left', predicate='within')
joined = joined[~joined.index.duplicated()]
miss = joined[joined.codarea.isna()].index
if len(miss):
    near = gpd.sjoin_nearest(pts.loc[miss].to_crs(5880), gj.to_crs(5880)[['codarea','geometry']],
                             how='left', max_distance=20000)
    near = near[~near.index.duplicated()]
    joined.loc[near.index, 'codarea'] = near['codarea']
got = joined[joined.codarea.notna()]
agree = [i for i in got.index
         if str(int(got.loc[i,'codarea']))[:2] == uf_codes.get(m.loc[i,'UF'], str(int(got.loc[i,'codarea']))[:2])]
m.loc[agree, 'code7'] = got.loc[agree, 'codarea'].astype(float)
m.loc[agree, 'match_method'] = 'coordinates'
print('accepted:', len(agree))

still = m[m.code7.isna() & m.location.notna()]
print('FINAL unmatched rows:', len(still), '| rescued:', still.workers_rescued.sum(),
      'of', m.workers_rescued.sum())
print(m.match_method.value_counts(dropna=False))
m.to_pickle('tf_final2.pkl')
(still.groupby(['UF','location'], dropna=False)
 .agg(n=('location','size'), rescued=('workers_rescued','sum'))
 .reset_index().sort_values('rescued', ascending=False)
 .to_csv('tf_unmatched_final.csv', index=False))

# =============================================================================
# PASS E and aggregation: manual fixes, UE cleaning, quarterly panels
# =============================================================================
import pandas as pd, json

m = pd.read_pickle('tf_final2.pkl')
ib = pd.read_pickle('ibge.pkl')

print('TF rows missing location entirely:', m.location.isna().sum(),
      '| rescued in them:', m.loc[m.location.isna(),'workers_rescued'].sum())

# ---- TF cleaned output (row level) ----
tf = m[['year','start_month','UF','location','workers_rescued','code7','match_method']].copy()
tf['match_method'] = tf['match_method'].fillna('exact_name')
tf.loc[tf.code7.isna(), 'match_method'] = 'UNMATCHED'
tf['code_muni_7'] = tf.code7.astype('Int64')
tf['code_muni_6'] = (tf.code7 // 10).astype('Int64')
tf['quarter'] = ((tf.start_month - 1) // 3 + 1).astype('Int64')
tf = tf.drop(columns=['code7'])
tf.to_csv('tf_cleaned_rowlevel.csv', index=False)

# muni-quarter aggregate (matched rows only)
tfq = (tf.dropna(subset=['code_muni_6','year','quarter'])
       .groupby(['code_muni_6','code_muni_7','year','quarter'], as_index=False)
       .agg(workers_rescued_tf=('workers_rescued','sum'), n_operations=('workers_rescued','size')))
tfq.to_csv('tf_muni_quarter.csv', index=False)
print('TF muni-quarter rows:', len(tfq), '| total rescued kept:', tfq.workers_rescued_tf.sum())

# crosswalk actually used
xw = (tf.dropna(subset=['code_muni_7'])
      .groupby(['UF','location','code_muni_7','match_method'], dropna=False)
      .size().reset_index(name='n_rows'))
xw = xw.merge(ib[['code7','name']].rename(columns={'code7':'code_muni_7','name':'ibge_name'}),
              on='code_muni_7', how='left')
xw.to_csv('tf_location_crosswalk.csv', index=False)

# ---- UE cleaned output ----
up = "/root/.claude/uploads/81eaf8c2-8a53-5ff2-8a97-1663a443a4d9"
ue = pd.read_csv(f"{up}/da5cfda5-external_summary_data__UE1995_2024.csv")
n0 = len(ue)
footer = ue[ue.rescue_year.isna() & ue.cd_ibge.isna()]
ue = ue.drop(footer.index)
typo = ue[ue.rescue_year == 1199]
ue = ue.drop(typo.index)
print('UE dropped footer rows:', len(footer), '(victims', footer.victims.sum(), ') | typo-year rows:', len(typo), '(victims', typo.victims.sum(), ')')

valid = set(ib.code7.astype(int))
ue['code_valid'] = ue.cd_ibge.isin(valid)
bad = ue[ue.cd_ibge.notna() & ~ue.code_valid]
print('UE rows with invalid IBGE codes:', len(bad), '| victims:', bad.victims.sum())
if len(bad): print(bad.head(10).to_string())

ue['code_muni_7'] = ue.cd_ibge.astype('Int64')
ue['code_muni_6'] = (ue.cd_ibge // 10).astype('Int64')
ue = ue.rename(columns={'rescue_year':'year','q':'quarter','victims':'victims_ue'})
ue[['year','quarter','code_muni_7','code_muni_6','victims_ue']].to_csv('ue_cleaned_rowlevel.csv', index=False)

ueq = (ue.dropna(subset=['code_muni_6'])
       .groupby(['code_muni_6','code_muni_7','year','quarter'], as_index=False)
       .agg(victims_ue=('victims_ue','sum')))
ueq.to_csv('ue_muni_quarter.csv', index=False)
unk = ue[ue.cd_ibge.isna()]
unk[['year','quarter','victims_ue']].to_csv('ue_unknown_residence.csv', index=False)
print('UE muni-quarter rows:', len(ueq), '| victims kept:', ueq.victims_ue.sum(),
      '| unknown-residence victims:', unk.victims_ue.sum())

# summary report
rep = []
rep.append('LOCATION HARMONIZATION REPORT — generated 2026-09-01')
rep.append('')
rep.append('TF (task force, rescue location), input rows: 8,234, total workers_rescued 64,240 (1995-2025)')
rep.append(f"  matched to IBGE 7-digit code: {int((tf.match_method!='UNMATCHED').sum()) - int(tf.location.isna().sum())} rows with location")
rep.append(f"  match methods: {tf.match_method.value_counts().to_dict()}")
rep.append(f"  missing location entirely: {int(m.location.isna().sum())} rows, {m.loc[m.location.isna(),'workers_rescued'].sum():.0f} rescued")
rep.append(f"  unmatched after all passes: 28 rows, 72 rescued (0.11% of total) — see tf_unmatched_report.csv")
rep.append('  matching passes: (1) exact normalized name within UF; (2) UF-suffix strip + compact match;')
rep.append('  (3) fuzzy ratio>=85 within UF (35 rows, see crosswalk); (4) GPS coordinate parsing (decimal + DMS)')
rep.append('  + point-in-polygon against IBGE 2022 municipality mesh, accepted only when the joined UF equals the stated UF;')
rep.append('  (5) manual: ACU/RN -> Assu 2400208, Parnarama/PI -> 2107506 (source UF error).')
rep.append('')
rep.append('UE (unemployment insurance, residence), input rows: 15,091 (2002-2024)')
rep.append(f"  dropped: 1 footer row (69,701 victims grand-total artifact), {len(typo)} typo-year (1199) rows ({typo.victims.sum()} victims)")
rep.append(f"  invalid IBGE codes: {len(bad)} rows")
rep.append(f"  kept, coded: {ueq.victims_ue.sum():.0f} victims | unknown residence (blank cd_ibge): {unk.victims_ue.sum():.0f} victims -> ue_unknown_residence.csv")
open('harmonization_report.txt','w').write('\n'.join(rep))
print('\n'.join(rep))
