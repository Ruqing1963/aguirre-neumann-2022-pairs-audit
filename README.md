# A note on the existence of pairs of *r*-primitive and *k*-normal elements over finite fields

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24313540.svg)](https://doi.org/10.5281/zenodo.24313540)

Paper, verification program and data for an audit of

> J. J. R. Aguirre and V. G. L. Neumann, *Pairs of r-primitive and k-normal elements in finite fields*,
> [arXiv:2210.11504](https://arxiv.org/abs/2210.11504) (2022);
> Bull. Braz. Math. Soc. (N.S.) 54 (2023), Paper No. 24,
> [doi:10.1007/s00574-023-00341-z](https://doi.org/10.1007/s00574-023-00341-z).

Only the arXiv text (v1) was checked. No claim is made about the journal version.

## Findings

1. **Theorem 3.2 is false as stated.** For `q = 2`, `n = 19`, `r_i = 1`, `k_i = 0` and `F = x + 1` the
   hypothesis holds (`2^9.5 ≈ 724 ≥ 128`), but `Tr(α + 1) = Tr(α) + 1`, so `α` and `α + 1` are never both
   normal. The same happens for `F = x^2 + x`.
2. **Cause.** The proof never treats the character-sum branch with `η_{d2}` trivial and `ψ_2` nontrivial.
   For affine `F` these terms have size `q^n`, not `O(q^{n/2})`.
3. **Theorem 4.8 fails for affine `F`.** If `p | n`, or `p ∤ n` and `gcd(n, q − 1) = 1`, no `F = ax + b`
   maps a 2-normal element to a 1-normal one (module-structure obstruction). At `(q, n) = (5, 10)` an
   exhaustive search over `F_{5^10}` finds no pair for the 16 affine members of `Υ_5(2,1)`, and between
   3800 and 4780 admissible `α` for each of the other 560.
4. **Corrected statements.** Theorem 3.2′ adds a nondegeneracy condition (ND): `y_1 x + y_2 F(x)` is not
   of Artin–Schreier form `h^p − h + c`. Theorem 4.8′ restricts to non-affine `F`.
5. **Independent numerical gap.** Procedure TotalSieve prints the constant `6` where the sieve criterion
   needs `36`. Of the 709 pairs `(q, n)` with `5 ≤ q ≤ 400`, `8 ≤ n ≤ 20` in the range of Theorem 4.8,
   202 cannot be certified by the sieve, 144 of them even with the constant `6`. These pairs are not
   counterexamples, but they are not proved.
6. **Untreated case for `n = 8`.** The final step of Proposition 4.6 covers the 4090405 prime powers
   `q < 1.45·10^8` with `q ≡ 1 (mod 6)`, but not the 4089777 with `q ≡ 5 (mod 6)`.

## Contents

| Path | Description |
|---|---|
| `aguirre_neumann_primitive_normal_pairs_audit.tex` / `.pdf` | The note (4 pages) |
| `code/verify_aguirre_neumann_2022.py` | Verification program (Parts A–E, and the figure) |
| `data/results.json` | Output of the program |
| `figures/fig_aguirre_neumann_2022_audit.pdf` | Figure 1 of the note |

## Reproducing the computations

Requires Python 3 with `numpy`, `sympy` and `matplotlib`.

```
pip install numpy sympy matplotlib
cd code
python verify_aguirre_neumann_2022.py
```

This rewrites `data/results.json` and `figures/fig_aguirre_neumann_2022_audit.pdf`. Options:

- `--skip510` reuses Part B from an existing `data/results.json` instead of re-enumerating `F_{5^10}`.
- `--png preview.png` also writes a PNG preview of the figure.

| Part | What it checks | Where in the note |
|---|---|---|
| A | `q = 2`, `n = 19`: hypothesis of Theorem 3.2 versus exhaustive count, for `F = x + 1` and `F = x^2 + x` | Example 1 |
| B | `q = 5`, `n = 10`: exhaustive count for all 576 members of `Υ_5(2,1)` | Example 4, Figure 1(a) |
| C | Exact value of the omitted-branch character sum for `F = x + 1` in `F_{5^10}` (run inside Part B) | Section 3 |
| D | Exact optimum of the sieve criterion with constants 36 and 6, for `5 ≤ q ≤ 400`, `8 ≤ n ≤ 20` | Section 5, Figure 1(b) |
| E | Prime-power counts below `1.45·10^8` by residue mod 6 | Section 5 |

## Building the paper

```
pdflatex aguirre_neumann_primitive_normal_pairs_audit.tex
pdflatex aguirre_neumann_primitive_normal_pairs_audit.tex
```

## Citation

Ruqing Chen, *A note on the existence of pairs of r-primitive and k-normal elements over finite fields*,
November 2026. [doi:10.5281/zenodo.24313540](https://doi.org/10.5281/zenodo.24313540)

## License

Code (`code/`) is released under the MIT License. The paper, figure and data are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See [LICENSE](LICENSE).

## Contact

Ruqing Chen, GUT Geoservice Inc., Montreal, Canada — ruqing@hotmail.com
