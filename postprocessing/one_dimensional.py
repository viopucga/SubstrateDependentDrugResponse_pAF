"""Quantities retained from the 1D cable protocols."""
import numpy as np
from scipy.signal import find_peaks


def activation_time_max_dvdt(time_ms, vm_mv) -> float:
    time = np.asarray(time_ms, dtype=float)
    vm = np.asarray(vm_mv, dtype=float)
    if time.ndim != 1 or vm.ndim != 1 or time.size != vm.size or time.size < 2:
        raise ValueError("time and vm must be equal-length one-dimensional arrays")
    dt = np.diff(time)
    if np.any(dt <= 0):
        raise ValueError("time must increase strictly")
    return float(time[int(np.argmax(np.diff(vm) / dt))])


def distal_conduction_velocity(
    middle_activation_ms: float,
    distal_activation_ms: float,
    middle_position_cm: float = 2.49,
    distal_position_cm: float = 3.75,
) -> float:
    """Return the distal CV (CV2), which is the value used in the paper."""
    elapsed_s = (float(distal_activation_ms) - float(middle_activation_ms)) / 1000.0
    if elapsed_s <= 0:
        raise ValueError("distal activation must occur after middle activation")
    return (float(distal_position_cm) - float(middle_position_cm)) / elapsed_s


def s2_conducted(vm_mv, prominence_mv: float = 20.0) -> bool:
    """Historical distal-node criterion: more than one peak means S2 conduction."""
    peaks, _ = find_peaks(np.asarray(vm_mv, dtype=float), prominence=prominence_mv)
    return bool(peaks.size > 1)


def effective_refractory_period(coupling_intervals_ms, conducted) -> float:
    """Return the last blocked coupling interval immediately before conduction."""
    ci = np.asarray(coupling_intervals_ms, dtype=float)
    cond = np.asarray(conducted, dtype=bool)
    if ci.ndim != 1 or cond.ndim != 1 or ci.size != cond.size or ci.size == 0:
        raise ValueError("coupling intervals and conducted flags must be equal-length 1D arrays")
    order = np.argsort(ci)
    ci, cond = ci[order], cond[order]
    transitions = np.flatnonzero((~cond[:-1]) & cond[1:])
    if transitions.size == 0:
        return float("nan")
    return float(ci[int(transitions[-1])])