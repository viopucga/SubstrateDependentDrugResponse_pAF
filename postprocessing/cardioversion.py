"""Cardioversion classification from a pseudo-electrogram."""
import numpy as np
from scipy.signal import find_peaks


def classify_cardioversion(
    ecg_data,
    lead_column: int = 5,
    prominence: float = 10.0,
    minimum_height: float = -40.0,
    minimum_distance_samples: int = 40,
    last_active_sample: int = 9700,
):
    """Return termination and detected peaks using the historical criterion.

    ``lead_column=5`` refers to zero-based column 5 in a matrix containing time
    plus nine pseudo-ECG leads.
    """
    data = np.asarray(ecg_data, dtype=float)
    if data.ndim != 2 or data.shape[1] <= lead_column:
        raise ValueError("ecg_data does not contain the requested lead column")
    signal = data[:, lead_column]
    peaks, properties = find_peaks(
        signal, prominence=prominence, height=minimum_height,
        distance=minimum_distance_samples,
    )
    terminated = peaks.size == 0 or int(peaks[-1]) < int(last_active_sample)
    return {"cardioversion": bool(terminated), "peak_indices": peaks, "peak_properties": properties}