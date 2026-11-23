# Exact W and generalized-GHZ Bell-measurement optima

The note proves the exact value `sqrt(13)/3` for the two-source closed-network
Bell-state-measurement optimization in Chakrabarty and Kundu,
arXiv:2610.35838v1, Section III D, page 5, with the definitions of Proposition 1.
It includes all central parity indices, including identity indices and repeated
pairs. The white-noise threshold is exactly `3/sqrt(13)`; the single-central-node
restriction has value `2*sqrt(2)/3`.
The neighboring generalized-GHZ conjecture in Eq. (15) has exact value
`sqrt(1+sin(2*alpha)^2)`. The note extends Chakrabarty–Kundu (2026) and
recovers the single-node optimization from Horodecki et al. (1995), rather
than presenting that consequence as a new optimization theorem.

For independent source visibilities `v1,v2`, the W optimum is
`sqrt(13*v1*v2)/3`, with violation exactly when `v1*v2 > 9/13`.
For standard-frame W swapping, the jointly heralded mean concurrence is
`2/3`, as in Chakrabarty–Kundu, Section IV; retaining just one central
outcome gives mean concurrence `2/9`.

## Files

- `note.tex`: self-contained LaTeX source with an embedded bibliography.
- `note.pdf`: readable typeset note.
- `metadata.json`: bibliographic metadata and scientific coverage.
- `verify.py`: exact Pauli, polynomial, and finite rational-frame checks.
- `verification.json`: deterministic output of those checks.
- `verify_siblings.py`: exact generalized-GHZ and W swapping checks.
- `sibling-verification.json`: their deterministic output.
- `coverage.json`: exact source locators, statements and proof coverage.
- `requirements.txt`: exact symbolic-algebra dependency for the checker.
- `LICENSE`: CC BY 4.0 license notice.

## Build

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively, use a standard LaTeX installation with the listed packages:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is embedded, so no BibTeX run is necessary.

## Exact checks

Python 3.11 or later and SymPy 1.14.0 are sufficient:

```sh
python3 -m pip install -r requirements.txt
python3 verify.py verification.json
python3 verify_siblings.py sibling-verification.json
```

Expected: `status` is `pass`; all 64 Pauli tensor entries and two polynomial
residuals agree; all 36,864 ordered parity pairs over 144 rational-frame pairs
satisfy the sharp bound; the standard frame attains `13/9`.
The all-frame result is proved in the note, not inferred from a finite sample.
The sibling checker verifies 64 symbolic generalized-GHZ Pauli entries,
three polynomial identities and 32,768 exact parameter/frame/parity cases,
including the product and maximally entangled endpoints. It also verifies
all sixteen W Bell outcomes and the spin-flip eigenvalues of all four
one-node conditional states. The generalized-GHZ theorem has an all-parameter
written proof in the note.

For each finite pair, the checker uses Eq. (11) of the posing paper:
`Lambda(u,v) = (s + sqrt(s*s - 4*g*g))/2`, where
`s = dot(u,u) + dot(v,v)` and `g = dot(u,v)`.
With `B = 13/9`, the comparison is made without square-root rounding by checking
`2*B >= s` and `B*B - B*s + g*g >= 0` in rational arithmetic.

## Formal statement and scientific scope

The pinned Lean source is
https://github.com/the-omega-institute/trureturing/blob/a252cc46f0e145418a978217b8abd41a73fd820d/D5/S3/Quantum/MeasurementGeometry/WStateBellMeasurementOptimum.lean.
Its `D5.S3.Quantum.MeasurementGeometry.WStateBellMeasurementOptimum.result`
has type `claim`, the exact supremum for the trace-defined W state, all
orthonormal frames, and all four central indices. Its `pair_upper_bound`,
`attainment` and `squared_optimum` establish the upper bound, standard-frame
witness and squared value `13/9`. Noise, the single-node restriction,
generalized-GHZ exactness and swapping probabilities have the written proofs
in the note. The larger central projective-measurement class for W and the
open-network W optimum remain separate questions. For maximally entangled
GHZ sources, product-projective optimality already follows from
Chakrabarty–Kundu, Theorem 4(i) at n=2, with Eq. (14).

## Authors and license

Haobo Ma (ChronoAI Pte Ltd, Singapore), then Wenlin Zhang (National University
of Singapore; ChronoAI Pte Ltd, Singapore). Copyright 2026 Haobo Ma and Wenlin Zhang.
Licensed under Creative Commons Attribution 4.0 International (CC BY 4.0):
https://creativecommons.org/licenses/by/4.0/.
