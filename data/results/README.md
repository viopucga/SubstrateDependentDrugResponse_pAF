# Processed simulation results

The CSV files contain processed biomarkers and outcomes, without membrane-voltage fields, time series, videos or solver output.

## `baseline_results.csv`

One row represents one population, profile and potential S2 rotor position (`i`, `m` or `f`). `s2_coupling_interval_ms` stores the actual S2 coupling interval. `rotor_present=0` and `NaN` biomarker values mean that no rotor existed at that position or that the value was not applicable.

APD90 is the BCL 1000 calibration value selected from the candidate-population biomarkers by population and profile ID. The 1D outputs retained are ERP and distal CV (`erp_1d_ms` and `cv_1d_cm_s`). Distal CV is calculated between 2.49 and 3.75 cm using the difference between activation times, defined by the maximum membrane-potential upstroke.

`df_baseline_hz` is the median dominant frequency across the nine pseudoEGMs for that individual rotor simulation. `area_baseline_cm2` is the corresponding rotor-tip trajectory area in cm². Both values are `NaN` when no sustained rotor was available at that position.

`df_baseline_interval_0p5hz` assigns each individual rotor to the baseline-DF interval used for Figure 4 (4.1-4.5, 4.5-5.0, 5.0-5.5, 5.5-6.0, 6.0-6.5 or 6.5-7.3 Hz). `df_baseline_subgroup` assigns the whole profile from its maximum available i/m/f baseline DF: `low` <= 4.9 Hz, `mid` > 4.9 and <= 6.5 Hz, and `high` > 6.5 Hz.

## `drug_response_results.csv`

One row represents one population, profile, rotor position, drug and concentration. `cardioversion=1` indicates termination and `cardioversion=0` indicates nontermination. `profile_responder=1` means that at least one available i/m/f rotor terminated for that profile, drug and concentration. It is `NaN` when no rotor was available.

`baseline_area_cm2` and `drug_area_cm2` contain rotor-tip trajectory areas. Post-drug DF and area are reported only for rotors sustained until the end of the simulation: `drug_df_hz` and `drug_area_cm2` are `NaN` when `cardioversion=1` or when the measurement is not applicable. Their baseline values are retained when a baseline rotor was available.

The table retains baseline and post-drug ERP, CV, DF and rotor-area values. Percentage changes can be calculated from these absolute baseline and drug values.

## Profile-level dominant frequency

Individual i/m/f values are retained. The public analyses use:

- the maximum i/m/f baseline DF only for the 4.1-7.3 Hz inclusion filter and the low/mid/high DF subgroups;
- the median i/m/f DF for histograms and drug-induced DF changes.

Missing values are omitted when calculating these profile summaries.

## `exclusive_responders.csv`

This table contains the 169 profiles that respond exclusively to one of the three Figure 6 drugs at the `x1` concentration: flecainide 1.5 µM, vernakalant 10 µM, or tertiapin-Q 0.1 µM. It is generated deterministically from `profile_responder` using `postprocessing.responders.exclusive_responders`; the three drug columns preserve the binary response pattern and `exclusive_drug` names the only effective drug.
