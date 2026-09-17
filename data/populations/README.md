# Model populations

`population_10000.csv` contains the 10,000 candidate ionic profiles. `profile_id` is the original one-based row number and is the stable identifier used throughout the repository. `PoM_candidate_parameters.mat` is the unmodified MATLAB source table.

The 13 dimensionless columns scale GCaL, GK1, GKr, GKs, GKur, Gto, GKACh, GNa, GNaK, GNCX, GSK, GK2P and GSAC. Candidate factors span approximately 0.5 to 2.0.

`calibrated_pAF_population.csv` and `calibrated_IREpAF_population.csv` contain the 853 and 282 final profiles used in the manuscript. Legacy simulation folders called these populations `nSR` and `pAF`, respectively. Public names consistently use `pAF` and `IREpAF`.

`calibration_biomarkers.csv` contains the archived 0D biomarkers at BCL 1000 and BCL 500 for every final profile. Column names carry the `_bcl1000` or `_bcl500` suffix. `passes_calibration` records the result of applying `bmkrs_paper.mat`, including the BCL 500 repolarization-failure and EAD exclusions.

`source_simulation_id` is retained for provenance. It is only unique within `source_batch`; `profile_id` is the public key. IREpAF combines an initial batch with an additional, more permissively calibrated batch.

`calibration_ranges.csv` is a human-readable conversion of the unmodified `bmkrs_paper.mat`. The wide APA, maximum-upstroke and plateau-amplitude ranges were intentionally nonrestrictive.

Candidate generation is disabled in the normal workflow because the exact candidate population is supplied. The optional generator requires explicit activation and uses the original Latin-hypercube approach with bounds 0.5 to 2.0.
