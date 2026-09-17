# Modified Courtemanche 0D model

The MATLAB implementation is based on the Courtemanche-Ramirez-Nattel human atrial model and includes the study modifications documented in the source header.

Public condition names are:

- `pAF`: called `nSR` in legacy folders;
- `IREpAF`: called `pAF` in legacy folders.

Persistent AF (`peAF`) branches have been removed.

## Run one profile

From MATLAB, with the repository root as the current directory:

```matlab
addpath('model/0d')
result = simulate_profile('pAF', 1);          % 50 beats, BCL 1000 ms
result = simulate_profile('IREpAF', 153);
```

The function loads the unmodified `data/populations/PoM_candidate_parameters.mat`. It returns Vm, intracellular calcium and last-beat biomarkers without writing simulation files.

## Optional candidate generation

Exact reproduction uses the supplied candidate population. `population/GeneratePopulation_1.m` is disabled by default. A new Latin-hypercube population requires `allow_generation=true`, an explicit random seed and uses bounds 0.5-2.0. Such a population will not reproduce the supplied candidates because the original random seed is unavailable.

## Calibration

`calibrate_population.m` accepts BCL 1000 and BCL 500 biomarker tables and applies the ranges in `data/populations/bmkrs_paper.mat`, followed by exclusion of repolarization failures and EADs at BCL 500.

The final 853/282 cohorts also reflect downstream selection and simulation availability; they are not simply the complete output of the 0D mask.

## Validation status

- Renaming legacy nSR to pAF is numerically exact at the equation level.
- A 50-beat pAF test reproduced archived biomarkers within 0.12 ms maximum absolute timing difference.
- The current IREpAF code reproduces the later simulation batch closely except for the integration-step-sensitive maximum-upstroke measure, which is outside the public 1D result scope.
- Both IREpAF tissue batches map to the later 581-profile calibration. All 282 final profile IDs pass `bmkrs_paper.mat` when evaluated with the biomarker simulations associated with that calibration.
