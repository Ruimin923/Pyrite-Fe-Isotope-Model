[README.md](https://github.com/user-attachments/files/34952228/README.md)
# Single-scenario 10-parameter Fe-isotope model

This release replaces the old environment-specific SL / YX forward-simulation and sensitivity-analysis design with **one generic reaction system and one common 10-dimensional parameter distribution**. It is a new numerical experiment: it must **not** reuse the old SL/YX `results.mat` files.

## Model equations

The model is a closed, three-reservoir (`O` = Fe(III), `L` = dissolved Fe(II), `P` = pyrite) isotope mass-balance model. It tracks separate `54Fe` and `56Fe` inventories.

- `JD = f O` (DIR)
- `JA = (1-f) O` (AIR)
- `JP = Psi L` (pyrite precipitation)
- `JX = kappaOX L` (Fe(II) reoxidation)
- `dO/dtau = -O + kappaOX L`
- `dL/dtau = O - (Psi+kappaOX)L`
- `dP/dtau = Psi L`

`H` has been absorbed into effective kinetic ratios and **is not independently sampled**. Physical rate definitions: `f=KD/(KD+KA H)`, `Psi=KP H/(KD+KA H)`, `kappaOX=KX/(KD+KA H)`, `tau=(KD+KA H)t`.

## Default common prior distribution

| Input order | Parameter | Bound(s) | Distribution |
|---|---|---|---|
| 1 | f | 0–1 | uniform |
| 2 | Psi | 0.1–10 | loguniform |
| 3 | beta | 0.00001–0.4 | loguniform |
| 4 | deltaO | -0.3–0.5 per mil | uniform |
| 5 | deltaL | -1–1 per mil | uniform |
| 6 | epsDIR | -3–-2 per mil | uniform |
| 7 | epsAIR | -2.5–-0.5 per mil | uniform |
| 8 | epsPY | -3–-0.5 OR +0.5–+3 per mil | 50/50 mixture, uniform within branch |
| 9 | epsOX | +0.4–+3 per mil | uniform |
| 10 | kappaOX | 0.01–1 | loguniform |

The combined `Psi`, `beta`, and `kappaOX` bounds are the **unions of the preceding SL and YX ranges**, with `beta` set to loguniform because the unified range spans several orders of magnitude. These are **new shared independent priors**; this is **not** an equal-weight mixture or fit to observed data. You can edit all bounds and `P.dist` entries in **`model_parameters.m`**.

## Running

Use MATLAB R2024b or a compatible version with Statistics and Machine Learning Toolbox and `ode15s`.

```matlab
cd('PATH_TO_EXTRACTED_FOLDER')
run_all('quick')   % short end-to-end smoke test; output_quick/
run_all('full')    % complete 7-scramble analysis; output/
run_all('reuse')   % postprocess only; requires compatible full output
```

The default `run_all()` is `quick` for safety. Use `run_all('full')` for publication results. A change to `model_parameters.m` requires rerunning `full`; the reuse mode verifies the saved complete prior structure and model version.

## Sobol design and outputs

- Independent inputs: 10; Sobol base dimension: 20 (A and B matrices).
- For each of 7 scrambled replicates: `A`, `B`, and 10 `AB_i` hybrid matrices (`12N` model evaluations).
- `NperRep = [1024, 1024, 1024, 1024, 1024, 1024, 2048]` in full mode.
- Single-scenario full design: **8192 A-draw trajectories** and **98,304 primary forward solutions**.
- `ode15s` event stop: `P=0.95` (total initial Fe normalized to one), `RelTol=1e-8`, `AbsTol=1e-12`.
- Growth targets: `Ptarget=[0.50, 0.80, 0.95]`, common grid `q=.05:.005:.95` (181 points).
- Outputs: `Y1=delta_new(.95)-delta_new(.05)`, `Y2=sum(abs(diff(delta_new)))-abs(Y1)`, `Dplus=1[Y1>.05]`, `Dminus=1[Y1<-.05]`.
- Jansen first-order S1 / total ST estimates with approximate 95% Student-t intervals across 6 equal-N independently scrambled replicates. The separate 2048-sample replicate checks convergence; no pairwise interaction simulations.
- Conditional Sobol summaries: **all priors**, **epsPY<0**, **epsPY>0**. The conditional branches restrict the sign; epsPY magnitude **continues to vary**.
- If an output has variance <= 1e-12, its Sobol estimates are NaN/undefined, not forcibly zero.

## Requested sensitivity figures

Saved in `output/figures/` (or `output_quick/figures/`) as PNG, PDF and SVG:

- **Fig 2** (`Fig02_sobol_sorted_95CI`): Y1 S1/ST with 95% CIs, 3 panels (all priors, epsPY<0, epsPY>0); 10 inputs per panel, sorted by ST.
- **Fig 2b** (`Fig02b_Y2_sobol_sorted_95CI`): matching Y2 sensitivities.
- **Fig 3** (`Fig03_sensitivity_heatmap`): 4 response metrics (Y1, Y2, Dplus, Dminus), total-effect index ST across 3 Ptarget values and 10 inputs.
- **Fig 4** (`Fig04_sobol_convergence`): S1 and ST convergence for the top 3 full-prior Y1 parameters; error bars for equal-N scrambles; largest-N replicate separately marked.
- **Fig 5** (`Fig05_direction_probabilities_epsPY_sign_conditioned`): **4 x 3** layout with **10 parameter panels and 2 legend/explanation panels**. The `epsPY` panel displays direct sign/magnitude effects across its negative and positive ranges. Each of the **other 9 parameters** displays four curves: increase/decrease, each conditioned separately on epsPY<0 and epsPY>0. Solid lines = negative branch; dashed lines = positive branch. The other priors are marginalized, including the **magnitude** of epsPY within its fixed-sign branch. Logarithmic inputs use log bins and the x axes display **the full prior bounds**.

Since this release targets generic sensitivity analysis rather than a SL/YX-specific case study, the former **Fig 1 and Supplementary Fig S1** (SL/YX representative cases and two-environment isotope ensembles) are not generated by this release.

## Output data layout

- `parameters_used.csv`: actual sampled bounds and distribution labels.
- `samples.csv`: 28 numeric columns: [rep, sample index, ten inputs, 12 outputs (`Y1` 3 targets, `Y2` 3 targets, `Dplus` 3 targets, `Dminus` 3 targets), 4 cumulative reaction fluxes]. The **epsPY sample is column 10**; **Y1 at Ptarget=.80 is column 14**.
- `sobol.csv`: [rep, N, target index, output index, parameter index, S1, ST].
- `sobol_ci_sorted.csv` and `sobol_ci_sorted_Y2.csv`: [branch index 0/1/2, ST rank, parameter index, S1, CI low, CI high, ST, CI low, CI high, number of usable scrambles].
- `scenario_sensitivity_heatmap.csv`: [target index, output index, parameter index, mean S1, mean ST].
- `sobol_convergence.csv`: [rep, N, target index, output index, parameter index, S1, ST].
- `directional_responses.csv`: [parameter index, epsPY sign group (0=direct,1=negative,2=positive), bin, direction (1=increase,2=decrease), lower, upper, center, probability, count].
- `probabilities.csv`: [target index, sign group 0/1/2, count, mean Y1, increase probability, decrease probability, near-zero probability, mean Y2].
- `ensemble_1.mat`, ... `ensemble_7.mat`: A/B/Y/U and instantaneous/cumulative isotope trajectories, for potential later custom plots.
- `VALIDATION.txt`: report from MATLAB validation after successful quick/full run.

## Validation limits

The delivered forward ODE equations were independently tested using Python/SciPy on 32 parameter sets, including conservation and a matrix-exponential check (`INDEPENDENT_CHECKS.txt`). **MATLAB/Octave was not available here**, so the scripts, actual Sobol output and publication figures require a MATLAB quick run followed by a full run on your computer. Sensitivity results must not be interpreted as observationally fitted parameters or constraints independent of the chosen priors.
