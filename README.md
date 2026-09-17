# Substrate-dependent drug response in paroxysmal atrial fibrillation

This repository accompanies the PLOS Computational Biology manuscript on substrate-dependent antiarrhythmic drug response in pAF and IREpAF populations.

It contains the modified Courtemanche 0D model, candidate and calibrated population parameters, static 1D/2D tissue geometry and settings, processed simulation biomarkers, and reusable postprocessing functions. ELVIRA, its binaries, simulation-launch code, raw voltage traces, pseudo-ECG time series, and manuscript figure scripts are outside this release.

## Contents

- `model/0d/`: MATLAB model and population simulation/calibration functions.
- `data/populations/`: candidate parameters, calibrated pAF/IREpAF cohorts, and calibration limits.
- `data/geometry/`: 1D and 2D VTK meshes and their tissue settings.
- `data/results/`: baseline, drug-response, and exclusive-responder CSV files.
- `postprocessing/`: reusable functions for action potentials, 1D measurements, dominant frequency, cardioversion, responder classification, and rotor analysis.

See the README in each directory for formats and conventions. Python dependencies are listed in `requirements.txt`.

## License and citation

The authors' original code is released under the GNU General Public License v3.0 only (see `LICENSE`). Please cite the accompanying article using `CITATION.cff`. Its PLOS Computational Biology citation is provisional until publication details and a DOI are available. The Courtemanche model source retains its original provenance and copyright notices; rights for third-party contributions should be confirmed before public release.