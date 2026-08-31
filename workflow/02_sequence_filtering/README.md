# 02 Sequence-based filtering

Seven ordered filters, defined in `configs/keap1_filters.json`. Positions are Nrf2 residue
numbers; the peptide C-terminus is anchored at Leu84, so for length *L* residue *r* sits at
column *L* − (84 − *r*).

| step | target | rule |
|---|---|---|
| 1 | E79, E82 | one of E, R, H, K, Y |
| 2 | T80 | not A, S or V |
| 3 | G81 | must be G |
| 4 | res77 | not N |
| 5 | res83 | not H |
| 6 | whole peptide | reject if it contains DG, DS, DT or DN |
| 7 | whole peptide | reject if it contains ETGE |

## Input and output

Takes a ProteinMPNN FASTA or a plain one-sequence-per-line file. Writes the surviving
sequences and, optionally, the count remaining after each step.
