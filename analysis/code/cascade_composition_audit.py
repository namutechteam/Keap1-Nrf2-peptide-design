#!/usr/bin/env python3
"""
Cascade composition audit for the Keap1-Nrf2 peptide design workflow.

Tracks the residue identity at the Nrf2 recognition positions (77-83) at every
checkpoint of the design cascade for which a complete sequence set was retained,
and quantifies how the ETGE-exclusion filter (filter 7) was satisfied.

Inputs  : analysis/data/sequences/{01_generated,02_filtered}/<L>mer/...
          Keap1_Nrf2_45peptide_benchmark.csv
Outputs : table_composition_by_stage.csv
          table_etge_break_position.csv
          table_position79_by_stage.csv

Usage   : python cascade_composition_audit.py <MPNN_EXPORT_DIR> <BENCHMARK_CSV> [OUTDIR]
"""
import sys, os, csv, collections

LENGTHS = [10, 11, 12, 13, 14, 15, 16]
POSITIONS = [77, 78, 79, 80, 81, 82, 83]
NATIVE = {77: 'D', 78: 'E', 79: 'E', 80: 'T', 81: 'G', 82: 'E', 83: 'F'}
MOTIF = (79, 80, 81, 82)


def col(residue, length):
    """0-based index of Nrf2 residue number within a peptide of this length."""
    return residue - 84 + length - 1


def read_fasta(path):
    """ProteinMPNN .fa: first record is the native template, remainder are designs."""
    seqs = [ln.strip() for ln in open(path) if ln.strip() and not ln.startswith('>')]
    return seqs[0], seqs[1:]


def read_list(path):
    return [ln.strip() for ln in open(path) if ln.strip()]


def composition(pairs, position):
    n = len(pairs)
    c = collections.Counter(s[col(position, L)] for s, L in pairs)
    return n, c


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


def main(export_dir, benchmark_csv, outdir='.'):
    os.makedirs(outdir, exist_ok=True)
    stages = collections.OrderedDict()

    generated, unique, f5, f6, f7 = [], [], [], [], []
    for L in LENGTHS:
        template, designs = read_fasta(_generated(export_dir, L))
        uniq = list(dict.fromkeys([template] + designs))     # template counted, as in Table 1
        generated += [(s, L) for s in designs]
        unique    += [(s, L) for s in uniq]
        f5        += [(s, L) for s in read_list(_filtered(export_dir, L, f'4IFL_{L}_0317.txt'))]
        f6        += [(s, L) for s in read_list(_filtered(export_dir, L, f'4IFL_{L}_0317_DGfil.txt'))]
        f7        += [(s, L) for s in read_list(_filtered(export_dir, L, f'4IFL_{L}_0317_DGfil_ETGEfil.txt'))]

    bench = [r for r in csv.DictReader(open(benchmark_csv)) if r['round'] == 'first-round']
    final28 = [(r['sequence'], int(r['length_aa'])) for r in bench]

    stages['Generated']                       = generated
    stages['Unique sequences']                = unique
    stages['Sequence filters 1-5']            = f5
    stages['Sequence filter 6 (aspartimide)'] = f6
    stages['Sequence filter 7 (ETGE excluded)'] = f7
    stages['Structure-based selection']       = final28

    # ---- Table A: native-residue retention by stage -------------------------
    with open(os.path.join(outdir, 'table_composition_by_stage.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['stage', 'n'] + [f'{NATIVE[p]}{p}_percent' for p in POSITIONS])
        for name, data in stages.items():
            row = [name, len(data)]
            for p in POSITIONS:
                n, c = composition(data, p)
                row.append(round(100 * c[NATIVE[p]] / n, 2))
            w.writerow(row)
            print(f'{name:38s} n={len(data):>7}  ' +
                  '  '.join(f'{NATIVE[p]}{p}:{100*composition(data,p)[1][NATIVE[p]]/len(data):6.2f}%'
                            for p in POSITIONS))

    # ---- Table B: how filter 7 survivors broke the ETGE pattern -------------
    counts = collections.Counter()
    for s, L in f7:
        dev = tuple(p for p in MOTIF if s[col(p, L)] != NATIVE[p])
        counts[dev] += 1
    with open(os.path.join(outdir, 'table_etge_break_position.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['deviating_position(s)', 'n', 'percent_of_5118'])
        for k, v in sorted(counts.items(), key=lambda x: -x[1]):
            lbl = '+'.join(f'{NATIVE[p]}{p}' for p in k) if k else 'none'
            w.writerow([lbl, v, round(100 * v / len(f7), 2)])

    # ---- Table C: full position-79 and -82 composition by stage -------------
    with open(os.path.join(outdir, 'table_position79_by_stage.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['stage', 'position', 'n'] + list('ACDEFGHIKLMNPQRSTVWY'))
        for name, data in stages.items():
            for p in (79, 80, 82):
                n, c = composition(data, p)
                w.writerow([name, p, n] + [round(100 * c[a] / n, 3) for a in 'ACDEFGHIKLMNPQRSTVWY'])

    print(f'\nETGE-containing sequences entering filter 7: '
          f'{sum(1 for s, L in f6 if all(s[col(p, L)] == NATIVE[p] for p in MOTIF))} / {len(f6)}')
    print(f'Outputs written to {os.path.abspath(outdir)}')


if __name__ == '__main__':
    main(*sys.argv[1:])
