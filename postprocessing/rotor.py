"""Rotor phase-singularity detection and trajectory-area calculation."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.fftpack import fft, ifft
from scipy.signal import butter, filtfilt


def load_quadrilateral_elements(path):
    """Read the four 1-based node columns used by ELVIRA and return zero-based IDs."""
    values = pd.read_csv(Path(path), sep=r"\s+", skiprows=2, header=None).values
    if values.shape[1] < 6:
        raise ValueError("element file must contain metadata followed by four node columns")
    return _validate_elements(values[:, 2:6] - 1)


def _validate_elements(elements):
    arr = np.asarray(elements)
    if arr.ndim != 2 or arr.shape[1] != 4:
        raise ValueError("elements must have shape (n_elements, 4)")
    if not np.issubdtype(arr.dtype, np.integer):
        if not np.all(arr == np.round(arr)):
            raise ValueError("element node IDs must be integers")
        arr = arr.astype(np.int64)
    else:
        arr = arr.astype(np.int64, copy=False)
    # Public array input is zero-based; files are converted by the loader.
    if np.any(arr < 0):
        raise ValueError("element node IDs cannot be negative")
    return arr


def _historical_otsu(image):
    """Otsu implementation used by the original code (fixed histogram range 0..1)."""
    counts, _ = np.histogram(np.asarray(image).ravel(), bins=256, range=(0, 1))
    total = np.asarray(image).size
    weighted_total = np.dot(np.arange(256), counts)
    bg_weight = bg_sum = 0
    best_variance = 0.0
    threshold = 0
    for i in range(256):
        bg_weight += counts[i]
        if bg_weight == 0:
            continue
        fg_weight = total - bg_weight
        if fg_weight == 0:
            break
        bg_sum += i * counts[i]
        bg_mean = bg_sum / bg_weight
        fg_mean = (weighted_total - bg_sum) / fg_weight
        variance = bg_weight * fg_weight * (bg_mean - fg_mean) ** 2
        if variance > best_variance:
            best_variance = variance
            threshold = i
    return threshold / 255.0


def rotor_area_from_density(tip_density, element_size_cm: float = 0.03):
    """Calculate the largest-component area from a final 2D tip-density map."""
    density = np.asarray(tip_density)
    if density.ndim != 2 or density.size == 0:
        raise ValueError("tip_density must be a non-empty two-dimensional array")
    threshold = _historical_otsu(density)
    mask = density > threshold
    labels, count = ndimage.label(mask)
    largest = np.zeros_like(mask, dtype=bool)
    if count:
        sizes = np.bincount(labels.ravel())[1:]
        largest = labels == (int(np.argmax(sizes)) + 1)
    return {
        "rotor_area_cm2": float(np.count_nonzero(largest) * element_size_cm**2),
        "largest_component": largest,
        "otsu_threshold": float(threshold),
    }

def _phase_matrix(vm, sampling_hz, node_chunk_size):
    """Filter and transform chunks of nodes while storing phase as float32."""
    b, a = butter(4, [3 / (sampling_hz / 2), 50 / (sampling_hz / 2)], btype="bandpass")
    phase = np.empty(vm.shape, dtype=np.float32)
    shift = abs(float(np.mean(vm, dtype=np.float64)))
    padlen = 3 * (max(len(a), len(b)) - 1)
    for first in range(0, vm.shape[1], node_chunk_size):
        last = min(first + node_chunk_size, vm.shape[1])
        block = np.asarray(vm[:, first:last], dtype=float) + shift
        filtered = filtfilt(b, a, block, axis=0, padtype="odd", padlen=padlen).astype(np.float32)
        spectrum = fft(filtered, axis=0)
        n_time = filtered.shape[0]
        spectrum[n_time // 2 + 1:, :] = 0
        spectrum[1:n_time // 2, :] *= 2
        phase[:, first:last] = np.angle(ifft(spectrum, axis=0)).astype(np.float32)
    return phase


def rotor_tip_and_area(
    vm,
    elements,
    start_sample=0,
    end_sample=None,
    sampling_hz=1000.0,
    phase_tolerance=0.1,
    element_size_cm=0.03,
    grid_shape=None,
    node_chunk_size=512,
    return_tip_events=False,
):
    """Detect rotor tips and calculate the largest accumulated trajectory area.

    Parameters use explicit zero-based indexing and an end-exclusive interval.
    ``vm`` must be time × nodes. ``elements`` must contain zero-based node IDs.
    The output density counts tip detections per element over the chosen interval.
    """
    vm = np.asarray(vm)
    elems = _validate_elements(elements)
    if vm.ndim != 2:
        raise ValueError("vm must have shape (time, nodes)")
    if elems.size == 0 or np.max(elems) >= vm.shape[1]:
        raise ValueError("element connectivity is incompatible with the vm node dimension")
    if start_sample < 0:
        raise ValueError("start_sample must be non-negative")
    end = vm.shape[0] if end_sample is None else int(end_sample)
    if end <= start_sample or end > vm.shape[0]:
        raise ValueError("the requested sample interval is invalid")
    if node_chunk_size < 1:
        raise ValueError("node_chunk_size must be positive")

    phase = _phase_matrix(vm, sampling_hz, node_chunk_size)
    density = np.zeros(elems.shape[0], dtype=np.uint32)
    events = [] if return_tip_events else None
    n1, n2, n3, n4 = (elems[:, i] for i in range(4))
    for time_i in range(int(start_sample), end):
        differences = np.column_stack((
            phase[time_i, n2] - phase[time_i, n1],
            phase[time_i, n3] - phase[time_i, n2],
            phase[time_i, n4] - phase[time_i, n3],
            phase[time_i, n1] - phase[time_i, n4],
        ))
        differences[differences > np.pi] -= 2 * np.pi
        differences[differences < -np.pi] += 2 * np.pi
        detected = np.flatnonzero(np.abs(np.abs(differences.sum(axis=1)) - 2 * np.pi) < phase_tolerance)
        density[detected] += 1
        if return_tip_events and detected.size:
            events.extend((time_i, int(element_i)) for element_i in detected)

    if grid_shape is None:
        side = int(round(np.sqrt(elems.shape[0])))
        if side * side != elems.shape[0]:
            raise ValueError("grid_shape is required when the element count is not a perfect square")
        grid_shape = (side, side)
    if int(np.prod(grid_shape)) != elems.shape[0]:
        raise ValueError("grid_shape does not match the number of elements")

    density_2d = np.flip(np.transpose(density.reshape(grid_shape)), axis=0)
    historical_view = np.flipud(np.fliplr(density_2d.T))
    area_result = rotor_area_from_density(historical_view, element_size_cm)
    return {
        "rotor_area_cm2": area_result["rotor_area_cm2"],
        "tip_density": historical_view,
        "largest_component": area_result["largest_component"],
        "otsu_threshold": area_result["otsu_threshold"],
        "tip_events": events,
        "sample_interval": (int(start_sample), end),
    }