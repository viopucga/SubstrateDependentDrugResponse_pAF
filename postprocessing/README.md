# Postprocessing

This package contains side-effect-free calculations extracted from the original simulation workflows. It does not launch ELVIRA or modify simulation directories.

- `action_potential.py`: last-beat 0D biomarkers, including APD90.
- `one_dimensional.py`: distal CV and ERP decisions.
- `dominant_frequency.py`: nine-lead pseudo-ECG DF and profile reductions.
- `cardioversion.py`: termination classification from pseudo-ECG activity.
- `responders.py`: Figure 6 exclusive-response classes.
- `rotor.py`: phase singularities, accumulated tip trajectory, and rotor area.

Inputs and numerical conventions are described in the function docstrings and the README files under `data/`. Install the dependencies with `pip install -r requirements.txt` from the repository root.

All array APIs use explicit zero-based indices. Rotor time intervals are end-exclusive. ELVIRA element files are read with `load_quadrilateral_elements`, which converts their one-based node identifiers to zero-based indices.