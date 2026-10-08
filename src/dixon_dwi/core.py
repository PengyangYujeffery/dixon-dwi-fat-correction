"""Kernel calibration and signal correction on aligned arrays."""
import numpy as np
from scipy.optimize import nnls


def shift_axis(volume, shift, axis=0):
    volume = np.asarray(volume)
    if not isinstance(axis, (int, np.integer)) or not -volume.ndim <= axis < volume.ndim:
        raise ValueError("invalid array axis")
    axis %= volume.ndim
    if not isinstance(shift, (int, np.integer)):
        raise ValueError("shift must be an integer")
    out = np.zeros_like(volume)
    if abs(shift) >= volume.shape[axis]:
        return out
    src, dst = [slice(None)] * volume.ndim, [slice(None)] * volume.ndim
    n = volume.shape[axis]
    if shift >= 0:
        src[axis], dst[axis] = slice(0, n - shift), slice(shift, n)
    else:
        src[axis], dst[axis] = slice(-shift, n), slice(0, n + shift)
    out[tuple(dst)] = volume[tuple(src)]
    return out


def near_masks(fat, distance, axis=0):
    plus, minus = np.zeros_like(fat, dtype=bool), np.zeros_like(fat, dtype=bool)
    for offset in range(1, distance + 1):
        plus |= shift_axis(fat, -offset, axis)
        minus |= shift_axis(fat, offset, axis)
    return plus, minus


def validate_kernel(shifts, weights):
    shifts, weights = np.asarray(shifts), np.asarray(weights, dtype=float)
    if shifts.ndim != 1 or weights.shape != shifts.shape or not len(shifts):
        raise ValueError("kernel shifts and weights must be equal-length vectors")
    if not np.all(np.isfinite(shifts)) or not np.all(shifts == np.round(shifts)):
        raise ValueError("kernel shifts must be finite integers")
    if len(np.unique(shifts)) != len(shifts):
        raise ValueError("kernel shifts must be unique")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("kernel weights must be finite and nonnegative")
    return shifts.astype(int), weights


def estimate_fat(fat, shifts, weights, scale=1.0, axis=0):
    """Positive shifts move signal toward increasing array indices."""
    shifts, weights = validate_kernel(shifts, weights)
    fat = np.asarray(fat, dtype=float)
    if not np.isfinite(scale) or scale <= 0 or not np.all(np.isfinite(fat)):
        raise ValueError("fat and scale must be finite; scale must be positive")
    out = np.zeros_like(fat)
    for offset, weight in zip(shifts, weights):
        if weight > 0:
            out += weight * shift_axis(fat, int(offset), axis)
    return scale * out


def subtract_fat(dwi, fat_signal):
    """Subtract the same residual-fat estimate from every volume, including b=0."""
    dwi, fat_signal = np.asarray(dwi, dtype=float), np.asarray(fat_signal, dtype=float)
    if dwi.shape[:-1] != fat_signal.shape:
        raise ValueError("DWI spatial shape must match the fat signal")
    if not np.all(np.isfinite(dwi)) or not np.all(np.isfinite(fat_signal)):
        raise ValueError("signals must be finite")
    return dwi - fat_signal[..., None]


def nonattenuating(b0, d1000, body, ff, slices, pe_axis=0):
    fat = body & (ff >= 0.80)
    plus, minus = near_masks(fat, 12, pe_axis)
    lean = body & (ff <= 0.12) & ~plus & ~minus & (b0 > 0)
    selected = np.zeros(body.shape, dtype=bool)
    selected[:, :, slices] = True
    lean &= selected
    if not np.any(lean):
        raise ValueError("no interior lean voxels for attenuation estimation")
    attenuation = float(np.median(d1000[lean] / b0[lean]))
    if not np.isfinite(attenuation) or not 0 <= attenuation < 1:
        raise ValueError("interior attenuation must lie in [0, 1)")
    return (d1000 - attenuation * b0) / (1 - attenuation)


def calibration_stats(b0, d1000, fat, body, ff, scale, slices=None, pe_axis=0,
                      shifts=None, min_voxels=500):
    b0, d1000, fat, ff = [np.asarray(v, dtype=float) for v in (b0, d1000, fat, ff)]
    body = np.asarray(body, dtype=bool)
    if b0.ndim != 3 or any(v.shape != b0.shape for v in (d1000, fat, body, ff)):
        raise ValueError("calibration arrays must have the same 3D shape")
    if pe_axis not in (0, 1):
        raise ValueError("PE axis must be in-plane (0 or 1)")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")
    if any(not np.all(np.isfinite(v)) for v in (b0, d1000, fat)) or not np.all(np.isfinite(ff[body])):
        raise ValueError("calibration signals and in-body fat fractions must be finite")
    if np.any((ff[body] < 0) | (ff[body] > 1)):
        raise ValueError("fat fractions must lie in [0, 1]")
    if min_voxels < 1:
        raise ValueError("min_voxels must be positive")
    slices = list(range(b0.shape[2])) if slices is None else list(slices)
    if not slices or len(set(slices)) != len(slices) or any(k < 0 or k >= b0.shape[2] for k in slices):
        raise ValueError("slices must be unique valid indices")
    shifts = np.arange(-40, 41) if shifts is None else np.asarray(shifts)
    shifts, _ = validate_kernel(shifts, np.zeros(len(shifts)))
    target = nonattenuating(b0, d1000, body, ff, slices, pe_axis)
    n = len(shifts) + 1
    stats = {"XtX": np.zeros((n, n)), "Xty": np.zeros(n), "n": 0}
    for k in slices:
        mask = body[:, :, k] & (b0[:, :, k] > 0)
        if mask.sum() < min_voxels:
            continue
        columns = [shift_axis(fat[:, :, k], int(s), pe_axis)[mask] for s in shifts]
        design = np.column_stack(columns + [np.ones(int(mask.sum()))])
        y = target[:, :, k][mask] / scale
        stats["XtX"] += design.T @ design
        stats["Xty"] += design.T @ y
        stats["n"] += len(y)
    if not stats["n"]:
        raise ValueError("no slices have enough calibration voxels")
    return stats



def solve_mixed(st, n_nonneg):
    XtX, Xty = (st['XtX'], st['Xty'])
    ncol = XtX.shape[0]
    nfree = ncol - n_nonneg
    T = np.zeros((ncol, n_nonneg + 2 * nfree))
    T[:n_nonneg, :n_nonneg] = np.eye(n_nonneg)
    T[n_nonneg:, n_nonneg:n_nonneg + nfree] = np.eye(nfree)
    T[n_nonneg:, n_nonneg + nfree:] = -np.eye(nfree)
    M = T.T @ XtX @ T
    lam = 0.001 * float(np.mean(np.diag(XtX)))
    L = np.linalg.cholesky(M + lam * np.eye(M.shape[0]))
    z, _ = nnls(L.T, np.linalg.solve(L, T.T @ Xty))
    return T @ z


def calibrate(cases, shifts=None, pe_axis=0, slices=None, min_voxels=500):
    """Pool calibration cases; keep evaluation cases separate."""
    shifts = np.arange(-40, 41) if shifts is None else np.asarray(shifts)
    shifts, _ = validate_kernel(shifts, np.zeros(len(shifts)))
    pooled = None
    for case in cases:
        stats = calibration_stats(**case, shifts=shifts, pe_axis=pe_axis,
                                  slices=slices, min_voxels=min_voxels)
        if pooled is None:
            pooled = stats
        else:
            for key in pooled:
                pooled[key] += stats[key]
    if pooled is None:
        raise ValueError("at least one calibration case is required")
    beta = solve_mixed(pooled, len(shifts))
    weights = beta[:len(shifts)]
    mass = float(weights.sum())
    centroid = float(shifts @ weights / mass) if mass > 0 else 0.0
    return {"shifts": shifts.tolist(), "weights": weights.tolist(),
            "mass": mass, "centroid_pixels": centroid, "pe_axis": pe_axis}
