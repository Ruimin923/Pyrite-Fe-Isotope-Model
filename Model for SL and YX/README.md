# SL/YX pyrite Fe-isotope model

This MATLAB source package models three closed Fe reservoirs (O, reactive Fe(III); L, dissolved Fe(II); P, pyrite) with four Fe fluxes. It does **not** simulate sulfur isotope profiles, spatial transport, or observational fitting. The model contains **10 independent uncertain inputs**; **H is not sampled**. `model_parameters.m` is the only file you need to edit for input bounds and distribution choices.

## Model and normalization

Given the underlying physical rate coefficients `KD`, `KA`, `KP`, `KX` and the environmental modifier H, define:

- `KR = KD + KA*H`, `tau = KR*t`
- `f = KD / KR`
- `Psi = KP*H / KR`
- `kappaOX = KX / KR`

`H` is absorbed into these *effective normalized coefficients* and is **not an independent model input**.

The resulting fluxes and total-Fe ODEs are:

```text
JD = f*O               (DIR: Fe(III) -> Fe(II))
JA = (1-f)*O           (AIR: Fe(III) -> Fe(II))
JP = Psi*L             (PY:  Fe(II) -> pyrite)
JX = kappaOX*L         (OX:  Fe(II) -> Fe(III))

dO/dtau = -O + kappaOX*L
dL/dtau =  O - (Psi+kappaOX)*L
dP/dtau =  Psi*L
```

Initial inventories: `O=1-beta`, `L=beta`, `P=0`. Total normalized Fe=1, represented approximately as 54Fe+56Fe. Each reaction transfers the two isotopes consistently with `Rproduct=Rsource*exp(epsilon/1000)` and a declared reference ratio `Rref=15.698`.

## Environment-specific priors

**Only three inputs differ between SL and YX.** All bounds and distributions are editable in `model_parameters.m`.

| Input | SL | YX | Distribution |
|---|---|---|---|
| f | [0,1] | [0,1] | uniform |
| **Psi** | **[1,10]** | **[0.1,1]** | loguniform |
| **beta** | **[1e-5,0.02]** | **[0.05,0.4]** | loguniform SL, uniform YX |
| deltaO | [-0.3,0.5] per mil | same | uniform |
| deltaL | [-1,1] per mil | same | uniform |
| epsDIR | [-3,-2] per mil | same | uniform |
| epsAIR | [-2.5,-0.5] per mil | same | uniform |
| epsPY | [-3,-0.5] or [0.5,3] per mil | same | 50:50 mixture of two uniform intervals |
| epsOX | [0.4,3] per mil | same | uniform |
| **kappaOX** | **[0.01,1]** | **[0.1,1]** | loguniform |

**`Psi` is a ratio of effective rate coefficients, not the actual instantaneous PY/DIR+AIR flux ratio:** `JP/(JD+JA)=Psi*L/O`. No H column should appear in parameter CSVs, plots or sensitivity tables in this version.

## Requirements

MATLAB R2024b and Statistics and Machine Learning Toolbox (`sobolset`, `scramble`, `net`, `prctile`, `tinv`). This package has been statically inspected and independently tested for the Fe/isotope equations outside MATLAB, but has **not** been run end-to-end in MATLAB within the packaging environment. The user's MATLAB run of `quick` is the final runtime check.

## Running

Extract in a writable directory, set MATLAB Current Folder to it, check no older functions shadow this release:

```matlab
which run_all -all
which fe_forward -all
which model_parameters -all
run_all('quick')
```

`quick` writes into `output_quick/` with three scrambled N=32 replicates and an independent N=64 replicate. **Quick results are not publication-grade**. Check `output_quick/VALIDATION.txt`. For final results:

```matlab
run_all('full')
```

`full`: six independent scrambled N=1024 replicates for approximate 95% Student-t uncertainty intervals, plus one independent N=2048 sample check. A 10-parameter Sobol design requires A, B and ten AB_i matrices (12 forward simulations per Sobol point); across SL and YX, the primary integration count is **196,608**, excluding tolerance/grid tests. Each environment has **8,192 A trajectories** (4,096 in each epsPY sign under the 50:50 mixture). Each trajectory yields Ptarget 0.50, 0.80 and 0.95 windows. Samples are not observed pyrite crystals.

To reprocess final results, without recalculating the primary ensemble:

```matlab
run_all('reuse')
```

`reuse` compares the complete prior struct and a kinetic model-version string; it refuses to reuse ensembles from other priors, especially older 11-input H-containing variants. Changing *any* prior requires `quick`/`full` again.

## Important outputs

- `output/results.mat`, `samples.csv`, `sobol.csv`, `sobol_convergence.csv`, `precision.csv`, `conservation.csv`
- `output/parameters_used.csv` records **the actual bounds and distributions** (do not edit this generated file)
- `output/representative_cases.csv` contains selected illustrative SL/YX trajectories, selected only from the **first A replicate**, to match the reported directional signs if available. They are not fits to observed amplitudes.
- `output/ensemble_pyrite_bands_targets.csv` shows **pointwise central 95% prior-predictive intervals** of instantaneous new and cumulative bulk pyrite delta56Fe, split by sign of epsPY and three Ptargets.
- `output/sobol_ci_sorted.csv` (Y1) and `sobol_ci_sorted_Y2.csv` (Y2) include S1/ST, approximate 95% scramble t intervals, ranked by ST.
- `output/VALIDATION.txt` is written only after all checks pass.

`Y1 = delta56Fe_new(q=.95) - delta56Fe_new(q=.05)`. `Y2 = total variation - abs(Y1)` measures nonmonotonicity, not geometric curvature. Dplus/Dminus use ±0.05 per mil thresholds. Index confidence intervals quantify **numerical scramble uncertainty**; prior-predictive bands describe **variation across input priors**, not measurement errors.

## Figures

| Name | Description |
|---|---|
| `Fig01_representative_trajectories` | SL/YX reservoir evolution and direction-compatible instantaneous/cumulative pyrite isotope examples, plus all-A isotope bands separated by epsPY sign at Ptarget=.80 |
| `Fig02_sobol_sorted_95CI` | Y1 S1/ST with scramble CIs, ordered by ST, for all/negative/positive epsPY |
| `Fig02b_Y2_sobol_sorted_95CI` | As Figure 2 but Y2 |
| `Fig03_sensitivity_heatmap` | ST across targets, environments and response metrics |
| `Fig04_sobol_convergence` | Variation with Sobol N |
| `Fig05_direction_probabilities_f_beta_epsPY_Psi_kappaOX` | **2 x 5** conditional probabilities vs f, beta, epsPY, Psi and kappaOX (H omitted) |
| `FigS01_Ptarget_isotope_envelopes` | Instantaneous/cumulative isotope bands for Ptarget=.50/.80/.95 |

The dedicated pairwise-interaction module and its heatmap are **absent**. Figure export: PNG, PDF and SVG; complex translucent-band plots use high-resolution PDF images to avoid MATLAB vector export warnings.

## Files

`run_all.m`, `run_project.m`, `fe_forward.m`, `analyze_results.m`, `make_figures.m`, `validate_outputs.m`, `model_parameters.m`. The initial sample columns have been **reindexed** for 10 inputs; old results cannot be used. `output/` folders are generated only by actually running MATLAB, not included here.

## Interpretation caveats

Results are prior-predictive, not posterior probabilities or experiments on natural mechanism prevalence. Rates/time are dimensionless. The O/L/P boxes are not explicit sediment-depth profiles. SL/YX differences are **assumed priors** on Psi, beta and kappaOX; they are not inferred from observed datasets. Because the normalization scale KR differs in principle between natural settings, comparison of dimensionless trajectories does not establish absolute physical reaction-time differences.


## Figure 5 revision

Figure 5 is now a 4 x 3 layout: SL/YX epsPY, f and Psi in the first two rows; SL/YX beta and kappaOX in rows three and four. For f, Psi, beta and kappaOX, probabilities are recalculated separately within the epsPY < 0 and epsPY > 0 subsets.
