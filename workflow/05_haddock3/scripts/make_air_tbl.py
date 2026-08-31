#!/usr/bin/env python3
"""Generate the HADDOCK3 AIR restraint file for a peptide of any length.

The restraints are defined once, in Nrf2 residue numbers, and translated to
peptide column indices for the requested length using the same C-terminal
anchor (Leu84) as the sequence filters in step 02:

    column = length - (84 - residue_number)

Usage:
    python make_air_tbl.py --length 15 --output 15_mer.tbl
    python make_air_tbl.py --self-test
"""

import argparse
import sys

ANCHOR_RESIDUE = 84            # Nrf2 Leu84, the fixed C-terminal residue
RECEPTOR_SEGID = "A"           # Keap1 Kelch domain
PEPTIDE_SEGID = "B"

# (active side, active id, passive side, [passive ids]) with peptide ids as Nrf2 numbers.
RESTRAINTS = [
    ("receptor", 380, "peptide", [82, 83]),
    ("receptor", 415, "peptide", [79]),
    ("receptor", 483, "peptide", [79]),
    ("peptide", 78, "receptor", [525, 530]),
    ("peptide", 79, "receptor", [415, 483, 508, 555]),
    ("peptide", 80, "receptor", [602]),
    ("peptide", 82, "receptor", [380, 363]),
]

DISTANCES = "2.0 2.0 0.0"


def column_of(residue_number, length, anchor=ANCHOR_RESIDUE):
    return length - (anchor - residue_number)


def _resid(side, number, length):
    """Receptor ids are Keap1 numbering; peptide ids are converted to column index."""
    return number if side == "receptor" else column_of(number, length)


def _segid(side):
    return RECEPTOR_SEGID if side == "receptor" else PEPTIDE_SEGID


def build_tbl(length):
    blocks = []
    for act_side, act_id, pas_side, pas_ids in RESTRAINTS:
        active = f"assign (resi {_resid(act_side, act_id, length)} and segid {_segid(act_side)})"
        passive = [f"       (resi {_resid(pas_side, p, length)} and segid {_segid(pas_side)})"
                   for p in pas_ids]
        body = "\n        or\n".join(passive)
        blocks.append(f"{active}\n(\n{body}\n) {DISTANCES}\n")
    return "\n".join(blocks)


def self_test():
    """Check the residue-to-column mapping and the generated restraint blocks."""
    assert column_of(79, 15) == 10 and column_of(79, 16) == 11
    assert column_of(84, 15) == 15                     # anchor is the last residue

    tbl15, tbl16 = build_tbl(15), build_tbl(16)

    # A380 pairs with peptide E82/res83: columns 13,14 in a 15mer and 14,15 in a 16mer
    assert "assign (resi 380 and segid A)" in tbl15
    assert "(resi 13 and segid B)\n        or\n       (resi 14 and segid B)" in tbl15
    assert "(resi 14 and segid B)\n        or\n       (resi 15 and segid B)" in tbl16

    # peptide E79 -> the four Keap1 arginine/serine partners
    assert "assign (resi 10 and segid B)" in tbl15
    assert "assign (resi 11 and segid B)" in tbl16
    for keap1_res in (415, 483, 508, 555):
        assert f"(resi {keap1_res} and segid A)" in tbl15

    assert tbl15.count("assign") == len(RESTRAINTS) == 7
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--length", type=int, help="peptide length in residues (10-16)")
    ap.add_argument("--output", help="output .tbl path (default: stdout)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return
    if not args.length:
        ap.error("--length is required")

    tbl = build_tbl(args.length)
    if args.output:
        with open(args.output, "w") as f:
            f.write(tbl)
        print(f"wrote {args.output}")
    else:
        print(tbl, end="")


if __name__ == "__main__":
    sys.exit(main())
