#!/usr/bin/env python3
"""Composition at the HADDOCK3 checkpoint, recovered for the 10-15-residue
classes only (115 of 140 sequences). The filter-7 and final-selection sets are
restricted to the same six length classes so the comparison is like for like.

Usage: python haddock_checkpoint_composition.py <MPNN_EXPORT_DIR> <HADDOCK_XLSX> <DATASET_CSV> [OUT]
"""
import csv, json, os, re, sys, collections
import openpyxl

AA = list('ACDEFGHIKLMNPQRSTVWY')
NATIVE = {77: 'D', 78: 'E', 79: 'E', 80: 'T', 81: 'G', 82: 'E', 83: 'F'}
PEP = re.compile(r'^[ACDEFGHIKLMNPQRSTVWY]{9,17}$')
LENGTHS = range(10, 16)


def col(res, L):
    return res - 84 + L - 1


def read_haddock_checkpoint(path, sheet='summary_05'):
    """Read the docking checkpoint from the released CSV or the original workbook."""
    if str(path).lower().endswith('.csv'):
        with open(path, newline='') as fh:
            rows = list(csv.DictReader(fh))
        return [(r['sequence'].strip(), int(r['length_aa'])) for r in rows
                if r.get('sequence', '').strip()]

    wb = openpyxl.load_workbook(path, data_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    out = []
    for j, head in enumerate(rows[0]):
        L = int(re.match(r'(\d+)mer', str(head)).group(1))
        for r in rows[1:]:
            v = r[j]
            if isinstance(v, str) and PEP.match(v.strip()) and len(v.strip()) == L:
                out.append((v.strip(), L))
    wb.close()
    return out


def _generated(export_dir, L):
    """Path to the generated FASTA, in either the released or the original layout."""
    for p in (os.path.join(export_dir, '01_generated', f'{L}mer', f'4IFL_{L}.fa'),
              os.path.join(export_dir, f'{L}mer', 'raw_sequences', f'4IFL_{L}.fa')):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f'generated sequences for {L}mer not found under {export_dir}')


def _filtered(export_dir, L, name):
    """Path to a filtered sequence list, in either the released or the original layout."""
    for p in (os.path.join(export_dir, '02_filtered', f'{L}mer', name),
              os.path.join(export_dir, f'{L}mer', 'filtered_sequences', name)):
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f'{name} not found under {export_dir}')


def main(export_dir, haddock_xlsx, dataset_csv, out='si_table_S4_haddock.json'):
    f7 = []
    for L in LENGTHS:
        p = _filtered(export_dir, L, f'4IFL_{L}_0317_DGfil_ETGEfil.txt')
        f7 += [(l.strip(), L) for l in open(p) if l.strip()]

    h140 = read_haddock_checkpoint(haddock_xlsx)
    sel = [(r['sequence'], int(r['length_aa']))
           for r in csv.DictReader(open(dataset_csv))
           if r['round'] == 'first-round' and int(r['length_aa']) <= 15]

    stages = [('Sequence filter 7 (10–15 aa)', f7),
              ('HADDOCK3 checkpoint (10–15 aa)', h140),
              ('Structure-based selection (10–15 aa)', sel)]

    rows = []
    for label, data in stages:
        n = len(data)
        for pos in (79, 80, 82):
            c = collections.Counter(s[col(pos, L)] for s, L in data)
            rows.append({'stage': label, 'position': pos, 'n': n,
                         **{a: round(100 * c[a] / n, 3) for a in AA}})
        print(f'{label:42s} n={n:>5}  ' +
              '  '.join(f'{NATIVE[p]}{p}:'
                        f'{100*collections.Counter(s[col(p,L)] for s,L in data)[NATIVE[p]]/n:5.1f}%'
                        for p in (79, 80, 82)))

    json.dump(rows, open(out, 'w'))
    print('->', os.path.abspath(out))


if __name__ == '__main__':
    main(*sys.argv[1:])
