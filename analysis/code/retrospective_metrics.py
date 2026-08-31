#!/usr/bin/env python3
"""Relate every computed prioritization metric to the measured inhibition of the
28 first-round designs. Reports (i) the rank of the single active design under
each metric and (ii) Spearman coefficients at both screening concentrations with
Benjamini-Hochberg correction. With one active compound these are descriptive."""
import csv, sys
from scipy import stats

rows = [r for r in csv.DictReader(open(sys.argv[1])) if r['round'] == 'first-round']
out  = sys.argv[2] if len(sys.argv) > 2 else 'table_metric_activity.csv'

METRICS = [
    ('chai_iptm_mean',                  'higher', 'Chai-1 ipTM'),
    ('chai_aggregate_mean',             'higher', 'Chai-1 aggregate score'),
    ('chai_receptor_peptide_iptm_mean', 'higher', 'Chai-1 receptor-peptide ipTM'),
    ('chai_ptm_mean',                   'higher', 'Chai-1 pTM'),
    ('haddock_score',                   'lower',  'HADDOCK score'),
    ('dockq',                           'higher', 'DockQ'),
    ('bsa_A2',                          'higher', 'Buried surface area'),
    ('evdw',                            'lower',  'E_vdw'),
    ('eelec',                           'lower',  'E_elec'),
    ('edesolv',                         'lower',  'E_desolv'),
    ('fnat',                            'higher', 'Fnat'),
    ('i_rmsd_A',                        'lower',  'i-RMSD'),
    ('l_rmsd_A',                        'lower',  'l-RMSD'),
    ('il_rmsd_A',                       'lower',  'il-RMSD'),
    ('rmsd_A',                          'lower',  'RMSD'),
]
ids  = [r['peptide_id'] for r in rows]
i11  = ids.index('Comp11')
inh1 = [float(r['inhibition_1uM_pct']) for r in rows]
inh10= [float(r['inhibition_10uM_pct']) for r in rows]

res = []
for key, direction, label in METRICS:
    v = [float(r[key]) for r in rows]
    hi = direction == 'higher'
    better = sum(1 for x in v if (x > v[i11] if hi else x < v[i11]))
    ties   = sum(1 for x in v if x == v[i11])
    rank   = f'{better+1}' if ties == 1 else f'{better+1}-{better+ties}'
    top    = max(range(len(v)), key=lambda i: v[i] if hi else -v[i])
    r1, p1  = stats.spearmanr(v, inh1)
    r10,p10 = stats.spearmanr(v, inh10)
    sign = 1 if hi else -1
    res.append(dict(metric=label, favorable=direction, rank_of_active=rank,
                    top_ranked=ids[top], top_ranked_inhibition_1uM=inh1[top],
                    rho_1uM=round(sign*r1, 3), p_1uM=round(p1, 3),
                    rho_10uM=round(sign*r10, 3), p_10uM=round(p10, 3)))

for tag in ('1uM', '10uM'):                       # Benjamini-Hochberg within each concentration
    ps = sorted(range(len(res)), key=lambda i: res[i][f'p_{tag}'])
    m  = len(res); prev = 1.0
    for k in range(m-1, -1, -1):
        i = ps[k]
        prev = min(prev, res[i][f'p_{tag}'] * m / (k+1))
        res[i][f'q_{tag}'] = round(min(prev, 1.0), 3)

with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0].keys())); w.writeheader(); w.writerows(res)

print(f"{'metric':30s}{'rank':>8}{'top-ranked (inh%)':>26}{'rho1':>8}{'q1':>7}{'rho10':>8}{'q10':>7}")
print('-'*94)
for r in res:
    print(f"{r['metric']:30s}{r['rank_of_active']:>8}"
          f"{r['top_ranked']+' ('+str(int(r['top_ranked_inhibition_1uM']))+'%)':>26}"
          f"{r['rho_1uM']:>8.3f}{r['q_1uM']:>7.3f}{r['rho_10uM']:>8.3f}{r['q_10uM']:>7.3f}")
print('-'*94)
print(f"Minimum q-value across all {len(res)} metrics x 2 concentrations: "
      f"{min(min(r['q_1uM'], r['q_10uM']) for r in res):.3f}  -> none significant after BH correction")
print(f"P(the single active ranks first of 28 under a non-informative metric) = 1/28 = {1/28:.4f}")
