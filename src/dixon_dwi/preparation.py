"""Array-level preparation helpers; image loading is left to the caller."""
import numpy as np
from scipy import ndimage


def spgr(alpha_deg, tr_ms, t1):
    a, e1 = (np.radians(alpha_deg), np.exp(-tr_ms / t1))
    return np.sin(a) * (1 - e1) / (1 - np.cos(a) * e1)


def shift2d_nan(img, si, sj):
    out = np.full(img.shape, np.nan)
    n0, n1 = img.shape
    s0, d0 = (slice(0, n0 - si), slice(si, n0)) if si >= 0 else (slice(-si, n0), slice(0, n0 + si))
    s1, d1 = (slice(0, n1 - sj), slice(sj, n1)) if sj >= 0 else (slice(-sj, n1), slice(0, n1 + sj))
    out[d0, d1] = img[s0, s1]
    return out


def resample(vol, A_src, A_dst, shape, offsets):
    inv = np.linalg.inv(A_src)
    out = np.zeros(shape)
    ii, jj = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]), indexing='ij')
    for k in range(shape[2]):
        acc = np.zeros(shape[:2])
        for oi in offsets:
            for oj in offsets:
                for ok in offsets:
                    pts = np.stack([ii.ravel() + oi, jj.ravel() + oj, np.full(ii.size, k + ok), np.ones(ii.size)])
                    vox = (inv @ (A_dst @ pts))[:3]
                    acc += ndimage.map_coordinates(vol, vox, order=1, mode='constant', cval=np.nan).reshape(shape[:2])
        out[:, :, k] = acc / len(offsets) ** 3
    return out
