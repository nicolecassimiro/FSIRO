# FSIRO — Fuzzy Surprise Interior-point Radiotherapy Optimization

A fuzzy fluence-map optimization method for external-beam radiotherapy, based on
the **surprise function**, and solved by a **custom primal-dual interior-point
solver**. FSIRO models each anatomical structure as a fuzzy number and minimizes
the total surprise of the delivered dose, comparing **triangular** and
**trapezoidal** membership functions. MATLAB, no toolboxes required.

> Companion code for the manuscript *"Fuzzy Radiotherapy Optimization Based on the
> Surprise Function: A Comparative Study of Triangular and Trapezoidal
> Memberships"* (under review).

---

## Method

Let $A$ be the dose-influence matrix and $x \ge 0$ the beamlet intensities, so
the dose is $d = A x$. Each structure is described by a fuzzy number whose
membership $\mu(d) \in [0,1]$ expresses how acceptable a dose is (1 = ideal). The
**surprise function** of Neumaier,

$$S(\mu) = \left(\tfrac{1}{\mu} - 1\right)^2,$$

penalizes deviation from the ideal, growing without bound as $\mu \to 0$. FSIRO
solves

$$\min_{x}\ \sum_i S\!\big(\mu_i(A x)\big)
\quad \text{s.t.}\quad A x = \omega,\ \ b_1 \le \omega \le b_E,\ \ 0 \le x \le U,$$

where the fuzzy **support** $[b_1, b_E]$ is enforced as a hard constraint
(so the dose stays within the admissible range) and $\mu$ is a **triangular**
(3-point) or **trapezoidal** (4-point) membership.

### What this code contributes

- **A custom primal-dual interior-point solver** for the fuzzy-surprise model
  (the full KKT system, search directions and step control are implemented here,
  rather than delegated to an off-the-shelf solver).
- To obtain robust convergence on this ill-conditioned problem, the solver
  employs a **Mehrotra predictor-corrector** strategy and solves the Newton step
  **matrix-free** by preconditioned conjugate gradients (never forming the dense
  normal matrix).
- A **$C^1$-smoothed surprise function**: $S$ is regularized below a threshold
  $\mu_{\min}$ by a quadratic extension that matches value and first derivative,
  removing the pole and the kink that otherwise prevent convergence. The
  ill-conditioning of the raw surprise, and this fix, are discussed in the paper.
- A **protocol-driven design**: structures are located by name and the fuzzy
  numbers are built from a small protocol file, so no per-patient index editing
  is needed. New treatment sites are added by writing a protocol, not by touching
  the solver.

The method is evaluated against the reference plans of the public
[**TROTS**](https://sebastiaanbreedveld.nl/trots/) benchmark (head-and-neck and
liver cases).

---

## Repository structure

```
FSIRO.m                 Main solver (interior-point, both membership shapes)
dvh_metrics.m           Generic DVH metrics for one structure (D98, V_x, HI, ...)
run_all_site.m          Shared batch engine used by the run_all_* wrappers

protocols/
  protocol_HN.m         Head-and-Neck protocol (structures, goals, supports, limits)
  protocol_liver.m      Liver SBRT protocol

evaluation/
  evaluate_plan_HN.m    Clinical (ICRU-83 / QUANTEC) plan evaluation, head-and-neck
  evaluate_plan_liver.m Clinical (QUANTEC) plan evaluation, liver

examples/
  example_HN.m          Minimal usage example
  run_all_HN.m          Batch: all HN patients x both shapes, DVH vs TROTS
  run_all_liver.m       Batch: all liver patients x both shapes, DVH vs TROTS
```

The TROTS patient `.mat` files are **not** redistributed here; download them from
the [TROTS website](https://sebastiaanbreedveld.nl/trots/) and place the
`Head-and-Neck/` and `Liver/` folders alongside the code.

---

## Requirements

MATLAB (R2018b or newer recommended). No additional toolboxes — the DVH
percentiles are implemented locally, so the Statistics Toolbox is not required.

---

## Quick start

Because the code is organized in subfolders, first add them to the path:

```matlab
addpath(genpath('.'))    % run once from the repository root
```

Then:

```matlab
% single plan
P   = protocol_HN();
out = FSIRO('Head-and-Neck/Head-and-Neck_01.mat', 'triangular', P, true);
evaluate_plan_HN(out.dose, out.sizes, out.Dp);

% compare triangular vs trapezoidal
example_HN

% liver
out = FSIRO('Liver/Liver_01.mat', 'triangular', protocol_liver(), true);
evaluate_plan_liver(out.dose, out.sizes, out.Dp);
```

`out` contains the beamlet vector `x`, the voxel dose `dose`, per-structure mean
doses, iteration count, run time, the stopping certificate `cp`, and the stopping
`criterion`.

### Batch runs

```matlab
run_all_liver('/path/to/data')   % ~minutes; validates the pipeline
run_all_HN('/path/to/data')      % longer; runs overnight
```

Each writes one `.mat` per plan under `resultados/<site>/planos/` and a tidy
`<site>_resumo.csv` with one row per (patient, structure, method), where `method`
is `triangular`, `trapezoidal`, or `TROTS` (reference). Batch runs are
**resumable**: plans already computed are not re-optimized.

---

## Adding a treatment site

Copy `protocols/protocol_liver.m`, change the structure name patterns, clinical
goals and supports, and pass it to `FSIRO`. The prescription and target reference
dose are read automatically from the `.mat` file. Note that each OAR **support**
must sit above the achievable per-voxel dose (structures adjacent to the target
unavoidably receive high dose); clinical dose limits are applied in the DVH
evaluation, not as fuzzy supports.

---

## How to cite

If you use FSIRO, please cite the manuscript (currently under review):

> N. C. Cassimiro and A. R. L. Oliveira, *"Fuzzy Radiotherapy Optimization Based
> on the Surprise Function: A Comparative Study of Triangular and Trapezoidal
> Memberships"*, under review, 2025.

A machine-readable citation is provided in [`CITATION.cff`](CITATION.cff)
("Cite this repository" on GitHub). This entry will be updated with the full
reference and DOI upon acceptance.

---

## Authors

- **Nicole C. Cassimiro** — University of Campinas (Unicamp)
- **Aurelio R. L. Oliveira** — University of Campinas (Unicamp)

## License

Released under the [MIT License](LICENSE).
