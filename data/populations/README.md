# Model populations

`population_10000.csv` contains the 10,000 candidate ionic profiles. `profile_id` is the original one-based row number and is the stable identifier used throughout the repository. `PoM_candidate_parameters.mat` is the unmodified MATLAB source table.

The 13 dimensionless columns scale GCaL, GK1, GKr, GKs, GKur, Gto, GKACh, GNa, GNaK, GNCX, GSK, GK2P and GSAC. Candidate factors span 0.5 to 2.0.

`calibrated_pAF_population.csv` and `calibrated_IREpAF_population.csv` contain the 853 and 282 final profiles used in the manuscript.

`calibration_biomarkers.csv` contains the archived 0D biomarkers at BCL 1000 and BCL 500 for every final profile. Column names carry the `_bcl1000` or `_bcl500` suffix. `passes_calibration` records the result of applying `bmkrs_paper.mat`, including the repolarization-failure and EAD exclusions.

`calibration_ranges.csv` contains the upper and lower limit values for each calibration biomarkers. The wide APA, maximum-upstroke and plateau-amplitude ranges were intentionally nonrestrictive.

Candidate generation is disabled in the normal workflow because the exact candidate population is supplied. The optional generator requires explicit activation and uses the original Latin-hypercube approach with bounds 0.5 to 2.0.
