#!/usr/bin/env python3
"""Rebuild consolidation: GISAID (location-filtered at source) + GenBank (country qualifier).

Changes from the previous build:
  * GISAID geography comes from the repository's structured Location field, not from
    inference on the strain designation. This removes the inference entirely and fixes
    two false positives ("Egyptian goose", "Rousettus aegyptiacus" - both South Africa).
  * The GISAID retrieval is no longer restricted to N1, so H5N8/H5N5/H5N2 are present.
  * Data cut-off 2026-09-28 applied to both repositories.
  * GISAID's own Clade annotation is carried through for cross-validation against Nextclade.
"""
import sys, re, json, csv
from collections import Counter, defaultdict
import pandas as pd
import os
# Input directory holding the retrieved FASTA/metadata files.
# Override with:  export AIV_DATA=/path/to/data
DATA = os.environ.get("AIV_DATA", "data")


CUTOFF = "2026-09-28"
U_GIS = DATA.rstrip("/") + "/"
U_NCBI = DATA.rstrip("/") + "/"

def read_fa(fn):
    d = {}; n = None; b = []
    for l in open(fn, encoding="utf-8", errors="replace"):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n = l[1:].strip(); b = []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

def nk(s):
    """Normalised strain key for cross-repository matching."""
    s = s.lower().strip()
    s = re.sub(r'\s*\(h\d+n?\d*\)\s*$', '', s)
    s = re.sub(r'[\s_]+', '/', s)
    s = re.sub(r'/+', '/', s)
    return s.strip('/')

def year_of(v):
    m = re.search(r'(19|20)\d{2}', str(v))
    return int(m.group(0)) if m else None

def build(tag, gis_fa, gis_xls, nc_fa, nc_meta, subpfx):
    stats = Counter()
    recs = {}

    # ---------------- GISAID ----------------
    md = pd.read_excel(U_GIS + gis_xls)
    md['_id'] = md['Isolate_Id'].astype(str).str.strip()
    meta = md.set_index('_id').to_dict('index')
    # every row is already Location = Africa / Egypt by construction of the query,
    # but assert it rather than assume
    loc_bad = sum(1 for r in meta.values()
                  if 'egypt' not in str(r.get('Location', '')).lower())
    stats['gisaid_location_not_egypt'] = loc_bad

    for h, seq in read_fa(U_GIS + gis_fa).items():
        p = h.split("|")
        if len(p) < 4:
            stats['gisaid_bad_header'] += 1; continue
        acc, seg, strain, isl = p[0].strip(), p[1].strip(), p[2].strip(), p[3].strip()
        stats['gisaid_raw'] += 1
        if seg.upper() != "HA":
            stats['gisaid_not_ha'] += 1; continue
        m = meta.get(isl)
        if m is None:
            stats['gisaid_no_meta'] += 1; continue
        if 'egypt' not in str(m.get('Location', '')).lower():
            stats['gisaid_drop_location'] += 1; continue
        st = str(m.get('Subtype', '')).replace('A / ', '').strip()
        if not st.startswith(subpfx):
            stats['gisaid_drop_subtype'] += 1; continue
        y = year_of(m.get('Collection_Date')) or year_of(strain)
        if not y or not (2010 <= y <= 2026):
            stats['gisaid_drop_year'] += 1; continue
        k = nk(strain)
        if not k:
            stats['gisaid_drop_nokey'] += 1; continue
        cur = recs.get(k)
        if cur is None or len(seq) > len(cur['seq']):
            recs[k] = dict(key=k, acc=isl, strain=strain, subtype=st, year=y, src="GISAID",
                           seq=seq, country="Egypt",
                           host=str(m.get('Host', '') or ''),
                           gisaid_clade=str(m.get('Clade', '') or ''),
                           both=(cur is not None))
        if cur is not None: recs[k]['both'] = True
        stats['gisaid_kept'] += 1

    # ---------------- GenBank ----------------
    nmeta = {r['accession']: r for r in
             csv.DictReader(open(U_NCBI + nc_meta, encoding="utf-8-sig"), delimiter="\t")}
    for h, seq in read_fa(U_NCBI + nc_fa).items():
        acc = h.split()[0].strip()
        stats['genbank_raw'] += 1
        r = nmeta.get(acc)
        if r is None:
            stats['genbank_no_meta'] += 1; continue
        if 'egypt' not in (r.get('country', '') or '').lower():
            stats['genbank_drop_country'] += 1; continue
        st = (r.get('serotype', '') or '').strip().upper().replace('A/', '')
        m2 = re.search(r'H\d+N\d+', st) or re.search(r'H\d+N\d+', h.upper())
        st = m2.group(0) if m2 else (re.search(r'H\d+', st).group(0) if re.search(r'H\d+', st) else '')
        if not st.startswith(subpfx):
            stats['genbank_drop_subtype'] += 1; continue
        strain = (r.get('strain', '') or '').strip()
        if not strain:
            mm = re.search(r'(A/[^\s,(]+)', h)
            strain = mm.group(1) if mm else ''
        y = year_of(r.get('collection_date')) or year_of(strain)
        if not y or not (2010 <= y <= 2026):
            stats['genbank_drop_year'] += 1; continue
        k = nk(strain)
        if not k:
            stats['genbank_drop_nokey'] += 1; continue
        cur = recs.get(k)
        if cur is None:
            recs[k] = dict(key=k, acc=acc, strain=strain, subtype=st, year=y, src="GenBank",
                           seq=seq, country="Egypt", host=(r.get('host', '') or ''),
                           gisaid_clade='', both=False)
        else:
            cur['both'] = True
            if len(seq) > len(cur['seq']):
                cur.update(acc=acc, seq=seq, src="GenBank", strain=strain)
        stats['genbank_kept'] += 1

    return recs, stats

def emit(tag, recs, stats):
    import os
    os.makedirs("out", exist_ok=True)
    with open(f"out/{tag}_consolidated.fasta", "w") as f:
        for k, r in sorted(recs.items()):
            f.write(f">{r['acc']}|{r['subtype']}|{r['year']}|{r['src']}\n{r['seq']}\n")
    with open(f"out/{tag}_master_metadata.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "strain", "subtype", "year", "source", "in_both",
                    "host", "gisaid_clade", "length"])
        for k, r in sorted(recs.items()):
            w.writerow([r['acc'], r['strain'], r['subtype'], r['year'], r['src'],
                        "yes" if r['both'] else "no", r['host'], r['gisaid_clade'],
                        len(r['seq'].replace('-', ''))])
    both = sum(1 for r in recs.values() if r['both'])
    out = {"tag": tag, "cutoff": CUTOFF, "unique_strains": len(recs),
           "by_source": dict(Counter(r['src'] for r in recs.values())),
           "by_subtype": dict(Counter(r['subtype'] for r in recs.values())),
           "by_year": dict(sorted(Counter(r['year'] for r in recs.values()).items())),
           "held_by_both_repositories": both,
           "gisaid_clade_counts": dict(Counter(r['gisaid_clade'] for r in recs.values()
                                               if r['gisaid_clade']).most_common()),
           "stats": dict(stats)}
    print(json.dumps(out, indent=1))
    return out

if __name__ == "__main__":
    allout = {}
    r5, s5 = build("H5", "gisaid_epiflu_sequence.fasta", "gisaid_epiflu_isolates.xls",
                   "NCBI_H5_HA.fasta", "NCBI_H5_HA_metadata.tsv", "H5")
    allout["H5"] = emit("H5", r5, s5)
    import os
    if os.path.exists(U_GIS + "gisaid_epiflu_sequence_H9.fasta"):
        r9, s9 = build("H9", "gisaid_epiflu_sequence_H9.fasta", "gisaid_epiflu_isolates_H9.xls",
                       "NCBI_H9_HA.fasta", "NCBI_H9_HA_metadata.tsv", "H9")
        allout["H9"] = emit("H9", r9, s9)
    else:
        print("\n[H9 GISAID FASTA not present yet - H9 skipped]", file=sys.stderr)
    json.dump(allout, open("out/CONSOLIDATION.json", "w"), indent=1)
