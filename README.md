

> **Logistic Gene Regulatory Networks: A Modelling Framework Beyond Hill Functions**
> Ismail Belgacem (2026). *Mathematical Biosciences* (submitted, MBS-D-26-00594).

The core contribution is a **recursive De Morgan translation** that turns any Boolean
gene-regulatory network into a smooth system of ordinary differential equations built
from *products of increasing and decreasing logistic functions*. Compared with the usual
Hill-function encoding, the logistic kernels are globally `C^∞`, analytically invertible,
carry a strictly positive basal expression rate, and decouple threshold from steepness.
This repository provides the translator (in Mathematica and R) together with every script
needed to regenerate the figures and numerical results in the paper.

---

## The translator

Each Boolean update rule is mapped, through the recursive correspondence

```
Φ(¬f)      = 1 − Φ(f)
Φ(f₁ ∧ f₂) = Φ(f₁) · Φ(f₂)
Φ(f₁ ∨ f₂) = 1 − (1 − Φ(f₁))·(1 − Φ(f₂))     (De Morgan)
Φ(xᵢ)      = f⁺(xᵢ, θᵢ, λᵢ)   (activation kernel)
Φ(¬xᵢ)     = f⁻(xᵢ, θᵢ, λᵢ)   (repression kernel)
```

to a product of logistic factors, and every gene obeys `xᵢ' = κᵢ·Φᵢ(x) − γᵢ·xᵢ`. The steepness
is Hill-matched by `λᵢ = n/θᵢ`. Rules are reduced to minimal disjunctive normal form before
translation (the soft map is neither idempotent nor distributive, so the canonical form makes
`Φ` a function of the Boolean *function* rather than of how it is written).

| File | Language | Role |
|---|---|---|
| `BooleanToLogisticODE.m`  | Mathematica package  | `BooleanToOdeSystem`, `logisticp`, `logisticm` (context `` BooleanToLogisticODE` ``) |
| `BooleanToLogisticODE.nb` | Mathematica notebook | notebook form of the same package |
| `BooleanToLogisticODE.R`  | R                    | standalone R implementation (sourced by the R figure scripts) |
| `verify_translator.R`     | R                    | checks over random inputs that the general translator reproduces the hand-coded right-hand sides exactly and recovers the Boolean values at the vertices |

---

## Reproducing the figures

Every figure-generating program is named after the figure it produces.

| Program(s) | Figure | Paper section | Output |
|---|---|---|---|
| `Figure2_two_gene_oscillator.R`                 | Fig. 2 | §3.1 Two-gene oscillator             | `Oscillateur_original.png` |
| `Figure3_toggle_switch.R`                       | Fig. 3 | §3.2 Genetic toggle switch           | `toggle_switch.png` |
| `Figure4_Traynard_cell_cycle.R` / `.nb`         | Fig. 4 | §5.4 Cell-cycle simulation           | `Traynard.*` |
| `Figure5_Traynard_robustness.R` / `.wl` / `.nb` | Fig. 5 | §5.5 Parameter robustness (LHS)      | `traynard_robustness*.png` |
| `Figure6_Verlingue_T2DM_geroconversion.R`/`.nb` | Fig. 6 | §5.6 Geroconversion (type-2 diabetes)| `verlingue_t2dm.png` |
| `Figure7_oscillator_formulation_comparison.R`   | Fig. 7 | §6.2 Hill / logistic / Samuilik       | `oscillator_comparison.png` |
| `Figure8_repression_function_comparison.R`      | Fig. 8 | §6.3 Repression-function comparison   | `Coparaison_3.png` |

`Verlingue_normal_variant.R` / `.nb` — the non-diabetic (healthy) variant of the Verlingue
network, provided for contrast with the geroconversion figure; it is not a figure in the paper.
(Figure 1 is a schematic drawn in LaTeX/TikZ and has no generating program.)

---

## Requirements

**R** (≥ 4.0) with packages `deSolve` (ODE integration) and `lhs` (Latin-hypercube sampling):

```r
install.packages(c("deSolve", "lhs"))
```

**Mathematica / Wolfram Language** (≥ 13.0) for the `.nb` and `.wl` files.

## Running

**R.** Run from the repository directory so the relative `source()` path resolves:

```bash
Rscript Figure2_two_gene_oscillator.R
Rscript Figure3_toggle_switch.R
# ... etc.
Rscript verify_translator.R          # translator self-check (prints PASS)
```
`Figure4/5/6` and `Verlingue_normal_variant` and `verify_translator` call
`source("BooleanToLogisticODE.R")`.

**Mathematica.** Keep `BooleanToLogisticODE.m` in the same directory as the notebook, then
evaluate top to bottom. The notebooks and `Figure5_Traynard_robustness.wl` load the package
with `` Needs["BooleanToLogisticODE`"] `` (the context name is fixed, so only the package file
must be present on the path).

---

## Parameters

The two-gene oscillator uses `λ=3, κ₁=3, γ₁=0.25, κ₂=4, γ₂=0.5, θ₁=4, θ₂=3, x₁(0)=x₂(0)=1`;
the toggle switch uses the symmetric `κ=10, γ=1, θ=5`. The Traynard cell-cycle and Verlingue
geroconversion parameters are set at the top of the corresponding scripts. The robustness sweep
draws `κᵢ ∼ U(50,100), γᵢ ∼ U(0.25,2), θᵢ ∼ U(10,20)` over a 33-dimensional Latin-hypercube box.
No experimental datasets are generated; all parameter values are drawn from the published
literature.

---

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Belgacem2026LogisticGRN,
  author  = {Belgacem, Ismail},
  title   = {Logistic Gene Regulatory Networks: A Modelling Framework Beyond Hill Functions},
  journal = {Mathematical Biosciences},
  year    = {2026},
  note    = {Manuscript MBS-D-26-00594}
}
```

© 2026 Ismail Belgacem.
