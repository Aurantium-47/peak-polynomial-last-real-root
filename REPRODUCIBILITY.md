# Reproduction scope

Use Python 3.10 or newer for the accompanying verifier, with SymPy as specified in `requirements.txt`. The separate coefficient checker and package-integrity checker use the standard library. Assertions must remain enabled: do not use Python's `-O` or `-OO` options.

## File integrity

```sh
python check_package.py
```

This checks all manifest-listed files, the two uncompressed canonical JSONL stream hashes, and the 4,180/4,784 row counts. It does not recompute the polynomial or root-count assertions in those rows. `SHA256SUMS.json` excludes itself; `PUBLIC_FILES.txt` lists the complete distribution including the manifest.

## Recorded limited checks

The actual limited reproduction used Python 3.8.3 and SymPy 1.6.1 and completed with exit code zero. This is a recorded environment, not a claim that the historical Python 3.13.5/SymPy 1.14.0 environment was recreated.

| Check | Recorded coverage |
|---|---|
| Separate standard-library coefficient implementation | 88 nonempty sets: 79 original-target sets and 9 singletons; 835 transform-identity points; 880 rational moment/sign points; 432 exact kernel points; 46,233 permutations |
| Accompanying verifier with maximum peak 10 and kernel maximum 10 | 88 sets; 582 moment/sign points; 864 kernel/two-step points; exact root and derivative checks |
| Comparison with the historical archive | 88 polynomial records and 864 kernel records matched; permutation tables matched; 88 monomial polynomials agreed with the independently reconstructed polynomials |

The commands corresponding to these checks are:

```sh
python verify_peak.py --max 10 --kernel-max 10 --output reproduced_small
python verify_independent_small.py
```

The second command creates `independent_results.json` in the current working directory. The distributed reference copy is `verification/independent_small_results.json`, so the default output does not overwrite it. Its mathematical body is unchanged from the executed implementation; only the module docstring was edited for distribution. Both source hashes are recorded in `verification/limited_reproduction.json`.

`verification/limited_run/` contains the small verifier's exact output. Packaging copied these records and checked their identities; it did not execute either mathematical program again. The small results do not replace the general proof.

## Historical larger records

`data/summary.json` reports Python 3.13.5 and SymPy 1.14.0 and the following earlier computations:

- Maximum peak 18: 4,180 nonempty admissible sets, comprising 4,163 original Conjecture 6.6 cases and 17 singletons.
- 27,663 specified moment/adjoining/sign checks.
- 4,784 kernel points, with kernel maximum 24 and 92 zero cases.
- Direct enumeration of 46,233 permutations of lengths 1 through 8.

The complete mathematical scan was **not rerun for this edition**. The stored summary, stream hashes, and record counts do not independently certify every earlier assertion or its process execution. A user wishing to regenerate that specification can run:

```sh
python verify_peak.py --max 18 --kernel-max 24 --output regenerated_large
```

No resource or runtime bound for that larger command is certified here. The optional `--input-zip` mode concerns a separate exploratory input archive that is not distributed and is unnecessary for this proof or the reproduction commands above.

## Exact-data conventions

Each gzip file stores canonical JSON Lines. Arrays are ascending coefficients except `newton`, which uses the basis binomial(x,j). Rationals are reduced `[numerator, denominator]` pairs with positive denominator. Integer roots include multiplicity; residual polynomials are primitive integer multiples with positive leading coefficient.

The verifier removes integer factors in [-3,m], separately checks them, and then counts remaining roots above a. The removal interval is not assumed to contain all roots. Transform endpoint factors at 0 and 1 are likewise removed before interior root counting. Gzip bytes may vary across compression versions while the uncompressed stream hash remains the same.

## Manuscript build

Using an existing LaTeX installation and the packages declared in `proof.tex`:

```sh
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error proof.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error proof.tex
```

Compilation verifies document generation; it is not a mathematical proof-assistant check.
