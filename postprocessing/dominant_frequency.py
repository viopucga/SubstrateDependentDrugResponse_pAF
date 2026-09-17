"""Dominant-frequency analysis of ELVIRA pseudo-electrograms."""
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks, periodogram


def dominant_frequency(egm, sampling_hz: float = 1000.0):
    """Return historical filtered and raw DF values and intermediate spectra."""
    x = np.asarray(egm, dtype=float)
    if x.ndim != 1 or x.size < 32 or not np.all(np.isfinite(x)):
        raise ValueError("egm must be a finite one-dimensional signal")
    b1, a1 = butter(4, [40 / (sampling_hz / 2), 250 / (sampling_hz / 2)], btype="band")
    rectified = np.abs(filtfilt(b1, a1, x))
    b2, a2 = butter(4, 20 / (sampling_hz / 2))
    filtered = filtfilt(b2, a2, rectified)
    nfft = 2 ** int(np.ceil(np.log2(filtered.size)))
    freq, filtered_psd = periodogram(filtered, fs=sampling_hz, window="hamming", nfft=nfft)
    _, raw_psd = periodogram(x, fs=sampling_hz, window="hamming", nfft=nfft)
    # Preserve the original half-open slice [first f>0.5, last f<15).
    first = int(np.flatnonzero(freq > 0.5)[0])
    last = int(np.flatnonzero(freq < 15)[-1])

    def strongest(psd):
        local, _ = find_peaks(psd[first:last])
        if local.size == 0:
            return 0.0
        return float(freq[first + local[int(np.argmax(psd[first:last][local]))]])

    return {
        "df_filtered_hz": strongest(filtered_psd), "df_raw_hz": strongest(raw_psd),
        "filtered_signal": filtered, "frequency_hz": freq,
        "filtered_psd": filtered_psd, "raw_psd": raw_psd,
    }


def dominant_frequency_nine_leads(ecg_data, sampling_hz: float = 1000.0):
    """Compute each lead DF and their median from time + nine lead columns."""
    data = np.asarray(ecg_data, dtype=float)
    if data.ndim != 2 or data.shape[1] != 10:
        raise ValueError("ecg_data must contain one time column and exactly nine pseudo-ECG leads")
    values = np.array([dominant_frequency(data[:, i], sampling_hz)["df_filtered_hz"] for i in range(1, 10)])
    return {"lead_df_hz": values, "median_df_hz": float(np.median(values))}


def profile_df(s2_df_hz, reduction: str = "median") -> float:
    """Reduce available S2_i/m/f values; median is the paper outcome."""
    values = np.asarray(s2_df_hz, dtype=float)
    if np.all(np.isnan(values)):
        return float("nan")
    if reduction == "median":
        return float(np.nanmedian(values))
    if reduction == "maximum":
        return float(np.nanmax(values))
    raise ValueError("reduction must be 'median' or 'maximum'")