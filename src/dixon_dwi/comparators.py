"""Physics-shaped and baseline-offset comparators."""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import nnls
from .models import DIRS

KS = np.arange(-40, 41)


def phys_shape(peaks, wfs, conv, sign, blur):
    pxppm = wfs * conv / 3.4
    k = np.zeros(len(KS))
    lo = int(KS[0])
    for ppm, amp in peaks:
        d = sign * ppm * pxppm
        i0 = int(np.floor(d))
        w = d - i0
        for idx, wt in ((i0, 1 - w), (i0 + 1, w)):
            if KS[0] <= idx <= KS[-1]:
                k[idx - lo] += amp * wt
    if blur:
        k = gaussian_filter1d(k, blur, mode='constant')
    return k / k.sum()


def fit_mass(st, shape):
    XtX, Xty, nk = (st['XtX'], st['Xty'], len(KS))
    A = np.array([[shape @ XtX[:nk, :nk] @ shape, shape @ XtX[:nk, nk]], [shape @ XtX[:nk, nk], XtX[nk, nk]]])
    b = np.array([shape @ Xty[:nk], Xty[nk]])
    c, icpt = np.linalg.solve(A, b)
    return (max(float(c), 0.0), float(icpt))


def noddi_offset(proto, y, v1):
    n = len(y)
    idx = np.empty(n, dtype=int)
    for s in range(0, n, 5000):
        idx[s:s + 5000] = np.argmax(np.abs(v1[s:s + 5000] @ DIRS.T), axis=1)
    out = np.full((n, 3), np.nan)
    eps = np.full(n, np.nan)
    for di in np.unique(idx):
        A, meta = proto.dictionary(di)
        A1 = np.column_stack([A, np.ones(A.shape[0])])
        for r in np.where(idx == di)[0]:
            x, _ = nnls(A1, y[r])
            xa, xi, xd = (x[:-2], x[-2], x[-1])
            s = xa.sum()
            if s > 0:
                out[r] = (xa @ meta[:, 0] / s, xa @ meta[:, 1] / s, xi / (s + xi))
            eps[r] = xd
    return (out, eps)


def fit_offset(proto, S0):
    S = S0.astype(np.float64)
    b0n = S[:, proto.b0].mean(axis=1)
    ok = b0n > 0
    S = np.maximum(S, 0.001 * np.abs(b0n)[:, None])
    y = S / np.where(ok, b0n, 1.0)[:, None]
    _, _, v1 = proto.dti(y)
    nd, eps = noddi_offset(proto, y, v1)
    denom = 1.0 - eps
    okc = ok & np.isfinite(eps) & (denom > 0)
    yc = (y - np.where(okc, eps, 0.0)[:, None]) / np.where(okc, denom, 1.0)[:, None]
    yc = np.maximum(yc, 0.001)
    fa, md, _ = proto.dti(yc)
    met = {'FA': np.where(okc, fa, np.nan), 'MD': np.where(okc, md, np.nan), 'Vintra': np.where(ok, nd[:, 0], np.nan), 'ODI': np.where(ok, nd[:, 1], np.nan), 'Vcsf': np.where(ok, nd[:, 2], np.nan)}
    return (met, np.where(ok, eps, np.nan))


def physics_kernel(water_fat_shift, reconstruction_scale, polarity, mass, blur=0.0):
    """Polarity and mass require calibration in the studied B3 implementation."""
    if polarity not in (-1, 1) or mass < 0 or water_fat_shift <= 0 or reconstruction_scale <= 0 or blur < 0:
        raise ValueError("invalid physics-kernel parameters")
    if not np.all(np.isfinite([water_fat_shift, reconstruction_scale, mass, blur])):
        raise ValueError("physics-kernel parameters must be finite")
    peaks = [(34.0 / 63.87, 8.83), (-32.1 / 63.87, 2.22)]
    shape = phys_shape(peaks, water_fat_shift, reconstruction_scale, polarity, blur)
    if not np.all(np.isfinite(shape)):
        raise ValueError("fat peaks fall outside kernel support")
    return KS.copy(), mass * shape
