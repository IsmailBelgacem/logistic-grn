# Logistic Gene Regulatory Networks

**Logistic Gene Regulatory Networks: A Modeling Framework Beyond Hill Functions**
Ismail Belgacem (2026). *Mathematical Biosciences* (submitted, MBS-D-26-00594).

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
translation (the soft map is neither idempotent nor distributive, so the canonical form makes `Φ` is a function of the Boolean *function* rather than of how it is written).

| File | Language | Role |
|---|---|---|
| `BooleanToLogisticODE.m`  | Mathematica package  | `BooleanToOdeSystem`, `logisticp`, `logisticm` (context `` BooleanToLogisticODE` ``) |
| `BooleanToLogisticODE.nb` | Mathematica notebook | notebook form of the same package |
| `BooleanToLogisticODE.R`  | R                    | standalone R implementation (`minimizeDNF`, `booleanToLogisticPhi`, `booleanToLogisticField`; sourced by every R figure script below) |

---

## Reproducing the figures

Every figure-generating program is named after the figure it produces in the current
manuscript numbering.

| Program(s) | Figure | Paper section | Output |
|---|---|---|---|
| `Figure2_two_gene_oscillator.R`                     | Fig. 2 | §3.1 Two-gene oscillator                                   | `Oscillateur_original.png` |
| `Figure3_toggle_switch.R`                           | Fig. 3 | §3.2 Genetic toggle switch                                 | `toggle_switch.png` |
| `Figure4_Traynard_cell_cycle.R` / `.nb`             | Fig. 4 | §5.4 Cell-cycle simulation                                 | `Traynard.png` |
| `Figure5_Traynard_robustness.nb`                    | Fig. 5 | §5.5 Parameter robustness (Latin-hypercube sweep)          | `traynard_robustness.png` |
| `Figure6_Verlingue_T2DM_geroconversion.R` / `.nb`   | Fig. 6 | §5.6 Geroconversion, type-2-diabetes variant                | `VERLINGUE_T2DM_VARIANT.png` |
| `Figure7_Verlingue_normal_variant.R` / `.nb`        | Fig. 7 | §5.7 Normal-physiology variant (insulin-resistance edge removed) | `VERLINGUE_NORMAL_VARIANT.png` |
| `Figure8_oscillator_formulation_comparison.R`       | Fig. 8 | §6.2 Hill / logistic / Samuilik comparison on the oscillator | `oscillator_comparison.png` |
| `Figure9_repression_function_comparison.R`          | Fig. 9 | §6.3 Repression-function comparison                         | `Coparaison_3.png` |



---

## Requirements

**R** (≥ 4.0) with the `deSolve` package (ODE integration):

```r
install.packages("deSolve")
```

**Mathematica / Wolfram Language** (≥ 13.0) for the `.nb` files.

## Running

**R.** Run from the repository directory so the relative `source()` path resolves:

```bash
Rscript Figure2_two_gene_oscillator.R
Rscript Figure3_toggle_switch.R
Rscript Figure4_Traynard_cell_cycle.R
Rscript Figure6_Verlingue_T2DM_geroconversion.R
Rscript Figure7_Verlingue_normal_variant.R
Rscript Figure8_oscillator_formulation_comparison.R
Rscript Figure9_repression_function_comparison.R
```
Each of these sources `BooleanToLogisticODE.R`; run them from this directory (or otherwise
make sure that file is on the working path) so the `source()` call resolves. Each script
also prints, alongside the plot, an independent console check of every Boolean fixed point
or dynamical claim it makes (e.g. `checkFixed`/`checkFixedPoint` calls), so a successful,
warning-free run is itself a correctness check.

**Mathematica.** Keep `BooleanToLogisticODE.m` in the same directory as the notebook, then
evaluate top to bottom. The notebooks load the package with
`` Needs["BooleanToLogisticODE`"] `` (the context name is fixed, so only the package file
must be present on the path). Note that a notebook's *cached* output cells (displayed plots,
generated ODE systems, echoed initial conditions) are only refreshed by re-evaluating the
notebook (Evaluation ▸ Evaluate Notebook); the source cells are the ground truth between
evaluations.

---

## Parameters

The two-gene oscillator uses `λ=3, κ₁=3, γ₁=0.25, κ₂=4, γ₂=0.5, θ₁=4, θ₂=3, x₁(0)=x₂(0)=1`;
the toggle switch uses the symmetric `κ=10, γ=1, θ=5`. The Traynard cell-cycle and both
Verlingue geroconversion variants (T2DM and normal-physiology) draw their kinetic
parameters and initial conditions from `U(50,100)` (production), `U(0.25,2)` (degradation),
and `U(10,20)` (threshold), respectively, then freeze the realisation actually used; each
is set at the top of the corresponding script and reported in the paper's parameter
tables. The robustness sweep draws `κᵢ ∼ U(50,100), γᵢ ∼ U(0.25,2), θᵢ ∼ U(10,20)` over a
33-dimensional Latin-hypercube box. No experimental datasets are generated; all network
topologies are drawn from the published literature.

---

## Citation

If you use this code, please cite the paper:

```bibtex
@article{Belgacem2026LogisticGRN,
  author  = {Belgacem, Ismail},
  title   = {Logistic Gene Regulatory Networks: A Modeling Framework Beyond Hill Functions},
  journal = {Mathematical Biosciences},
  year    = {2026},
  note    = {Manuscript MBS-D-26-00594}
}
```
and, if you wish to cite the software archive directly, the Zenodo record
(`https://doi.org/10.5281/zenodo.21321885`).

## Links

- Paper repository: <https://github.com/IsmailBelgacem/logistic-grn>
- Software archive (Zenodo): <https://zenodo.org/records/21321885>

© 2026 Ismail Belgacem.
