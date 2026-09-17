"""Pure functions for biomarkers measured from a paced 0D action potential."""
import numpy as np
from scipy.signal import find_peaks


def select_time_window(time, values, start_ms: float, end_ms: float):
    """Return samples in the inclusive interval ``[start_ms, end_ms]``."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if time.ndim != 1 or values.ndim != 1 or time.size != values.size:
        raise ValueError("time and values must be one-dimensional arrays of equal length")
    if time.size == 0 or np.any(np.diff(time) < 0):
        raise ValueError("time must be non-empty and monotonically increasing")
    mask = (time >= start_ms) & (time <= end_ms)
    if not np.any(mask):
        raise ValueError("the requested time window contains no samples")
    return time[mask], values[mask]


def action_potential_biomarkers(time_ms, vm_mv, peak_prominence_mv: float = 10.0):
    """Reproduce the biomarkers in the original ``get_AP_biomarkers`` routine.

    AP durations are measured from maximum dV/dt to the first sampled value at
    or below each repolarization threshold. No crossing interpolation is used,
    matching the historical implementation.
    """
    time = np.asarray(time_ms, dtype=float)
    vm = np.asarray(vm_mv, dtype=float)
    if time.ndim != 1 or vm.ndim != 1 or time.size != vm.size or time.size < 3:
        raise ValueError("time and vm must be equal-length 1D arrays with at least 3 samples")
    dt = np.diff(time)
    if np.any(dt <= 0) or not np.all(np.isfinite(time)) or not np.all(np.isfinite(vm)):
        raise ValueError("time must increase strictly and all samples must be finite")

    peaks, _ = find_peaks(vm, prominence=peak_prominence_mv)
    if peaks.size == 0:
        raise ValueError("no action potential with the requested prominence was found")
    peak_i = int(peaks[0])
    rmp = float(min(vm[0], vm[-1]))
    apo = float(vm[peak_i])
    apa = apo - rmp
    dvdt = np.diff(vm) / dt
    dvdt[~np.isfinite(dvdt)] = 0.0
    upstroke_i = int(np.argmax(dvdt))

    def apd(fraction: float):
        threshold = apo - apa * fraction
        candidates = np.flatnonzero(vm[peak_i:] <= threshold)
        idx = peak_i + int(candidates[0]) if candidates.size else vm.size - 1
        return float(time[idx] - time[upstroke_i]), idx

    apd95, _ = apd(0.95)
    apd90, apd90_i = apd(0.90)
    apd50, _ = apd(0.50)
    apd30, apd30_i = apd(0.30)
    apd20, apd20_i = apd(0.20)
    plateau_20 = float(np.mean(vm[apd20_i:apd30_i + 1]))
    plateau_time = time[peak_i] + 20.0
    plateau_candidates = np.flatnonzero(time >= plateau_time)
    plateau_i = int(plateau_candidates[0]) if plateau_candidates.size else vm.size - 1

    ead_peaks, _ = find_peaks(vm[peak_i:apd90_i + 1], prominence=2, width=5)
    dad_peaks, _ = find_peaks(dvdt[apd90_i:], prominence=0.1)
    return {
        "RMP": rmp, "APA": apa, "PlatA": float(vm[plateau_i] - rmp),
        "mx_dvdt": float(np.max(dvdt)), "APD95": apd95, "APD90": apd90,
        "APD50": apd50, "APD20": apd20,
        "rep_failure": int(vm[-1] > -40), "rep_EAD": int(ead_peaks.size > 1),
        "rep_DAD": int(dad_peaks.size > 0), "PLT_20": plateau_20,
        "triang": float(apd90 - apd30),
    }