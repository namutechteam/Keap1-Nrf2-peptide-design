#!/usr/bin/env python3
"""Merge Chai-1 confidence (mean over the five predicted models) into the
45-peptide dataset and write the machine-readable table released with the paper."""
import csv, statistics as st, sys, os

bench_csv, chai_all, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

by = {}
for r in csv.DictReader(open(chai_all)):
    by.setdefault(r['sequence'], []).append(r)

rows = list(csv.DictReader(open(bench_csv)))
extra = ['chai_iptm_mean', 'chai_ptm_mean', 'chai_receptor_peptide_iptm_mean',
         'chai_aggregate_mean', 'chai_n_models']
fields = list(rows[0].keys()) + extra

for r in rows:
    m = by.get(r['sequence'])
    if m:
        r['chai_iptm_mean']  = f"{st.mean(float(x['iptm']) for x in m):.4f}"
        r['chai_ptm_mean']   = f"{st.mean(float(x['ptm']) for x in m):.4f}"
        r['chai_receptor_peptide_iptm_mean'] = f"{st.mean(float(x['receptor_peptide_iptm']) for x in m):.4f}"
        r['chai_aggregate_mean'] = f"{st.mean(float(x['aggregate_score']) for x in m):.4f}"
        r['chai_n_models'] = len(m)
    else:
        for k in extra: r[k] = ''

with open(out_csv, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)

got = sum(1 for r in rows if r['chai_iptm_mean'])
print(f"{len(rows)} peptides written; Chai-1 confidence attached to {got}")
print("->", os.path.abspath(out_csv))
