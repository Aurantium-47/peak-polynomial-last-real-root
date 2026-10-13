# The last real root of a peak polynomial

Aurantium

Let S be a finite admissible peak set with at least two elements, m its largest element, and a its second-largest element. The manuscript gives a self-contained integral proof that

- p_S(x) < 0 for a < x < m;
- p_S(m) = 0;
- p_S(x) > 0 for x > m;
- p'_S(m) > 0.

Consequently m is the only real root above a, as asserted in Billey--Fahrbach--Talmage, arXiv:1410.8506v2, Conjecture 6.6. The argument imposes no bound on the number, positions, or gaps of the peaks. It also includes the auxiliary singleton case with a=0.

## Relation to earlier results

The focus is the complete interval (a,m) and the self-contained sign argument. The recurrence and finite-difference symmetry identity have earlier sources, including BFT Theorem 4.3. Li's version 3.0 coefficient theorem, DOI [10.5281/zenodo.22832569](https://doi.org/10.5281/zenodo.22832569), already implies right-side positivity, simplicity at m, negativity on [m-1,m), the positive symmetric transform properties, and the complete-monotonicity consequence. These are credited and reproved for completeness; they are not each presented as new discoveries.

The remaining noninteger portion of (a,m-1) requires the additional sign analysis developed in the manuscript. The literature comparison is bounded and does not certify global priority for the complete root-exclusion result. See REFERENCES.md and the manuscript's attribution paragraphs.

## Proof and computations

Read `proof.pdf` and its editable source `proof.tex`. Finite computations are diagnostic evidence, not premises of the all-parameter proof.

`data/` preserves historical records for 4,180 peak sets and 4,784 kernel points. The complete historical suite was not mathematically rerun for this edition. Hash and record-count checks confirm file consistency, not every historical execution.

The separately recorded limited reproduction covers maximum peak 10: 88 sets using the accompanying verifier, together with a separate exact coefficient implementation. The limited verifier's 88 polynomial and 864 kernel records match the corresponding archived records. The actual scope, commands, environments, and distinction from the larger archive are documented in REPRODUCIBILITY.md and `verification/limited_reproduction.json`. Packaging itself did not rerun the mathematical checks.

## Files

- `proof.tex`, `proof.pdf`: complete mathematical argument, attribution, and scope.
- `verify_peak.py`: exact finite verifier and data generator; its mathematical code is unchanged from the archived version.
- `verify_independent_small.py`: separate standard-library coefficient implementation; it writes `independent_results.json` in the current working directory.
- `check_package.py`: file hashes, exact-data stream hashes, and record counts only.
- `data/`: historical larger records and their summary.
- `verification/`: recorded limited results and execution-scope summary.
- `REPRODUCIBILITY.md`, `REFERENCES.md`, `RIGHTS.md`: reproduction, attribution, and distribution scope.
- `PUBLIC_FILES.txt`, `SHA256SUMS.json`: exact public file list and checksum manifest.

This work was developed with AI assistance. No external expert peer-review, proof-assistant certification, or exhaustive priority claim is made. The manuscript does not classify all other real roots or establish the separate complex-root conjectures.
