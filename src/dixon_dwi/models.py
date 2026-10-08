"""Diffusion model fitting used in the study."""
import numpy as np
from scipy.optimize import nnls

DPAR, DISO = 1.7, 3.0
VF_GRID = np.linspace(0.1, 0.99, 12)
OD_GRID = np.hstack(([0.03, 0.06], np.linspace(0.09, 0.99, 10)))


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


class Protocol:

    def __init__(self, bval, bvec):
        self.bv = bval / 1000.0
        g = bvec.T.astype(np.float64)
        n = np.linalg.norm(g, axis=1)
        g[n > 0] /= n[n > 0, None]
        g[bval <= 1] = 0.0
        self.g = g
        self.b0 = bval <= 1
        bv = self.bv
        self.B = np.column_stack([np.ones_like(bv), -bv * g[:, 0] ** 2, -bv * g[:, 1] ** 2, -bv * g[:, 2] ** 2, -2 * bv * g[:, 0] * g[:, 1], -2 * bv * g[:, 0] * g[:, 2], -2 * bv * g[:, 1] * g[:, 2]])
        self.Bpinv = np.linalg.pinv(self.B)
        self.estick = np.exp(-bv[:, None] * DPAR * (g @ SPHERE.T) ** 2)
        self.eiso = np.exp(-bv * DISO)
        self.cache = {}

    def dictionary(self, di):
        if di not in self.cache:
            mu = DIRS[di]
            gm2 = (self.g @ mu) ** 2
            c2 = (SPHERE @ mu) ** 2
            cols, meta = ([], [])
            for od in OD_GRID:
                kap = 1.0 / np.tan(od * np.pi / 2)
                w = np.exp(kap * (c2 - 1.0))
                w /= w.sum()
                e_ic = self.estick @ w
                tau = float(w @ c2)
                for vf in VF_GRID:
                    dperp = DPAR * (1 - vf)
                    da, dr = (dperp + (DPAR - dperp) * tau, dperp + (DPAR - dperp) * (1 - tau) / 2)
                    cols.append(vf * e_ic + (1 - vf) * np.exp(-self.bv * (dr + (da - dr) * gm2)))
                    meta.append((vf, od))
            cols.append(self.eiso)
            self.cache[di] = (np.column_stack(cols), np.array(meta))
        return self.cache[di]

    def dti(self, y):
        ly = np.log(np.maximum(y, 1e-06))
        beta = ly @ self.Bpinv.T
        for _ in range(2):
            w = np.exp(np.clip(2 * beta @ self.B.T, -50, 50))
            G = np.einsum('nm,mi,mj->nij', w, self.B, self.B)
            h = np.einsum('nm,mi,nm->ni', w, self.B, ly)
            beta = np.linalg.solve(G + 1e-12 * np.eye(7), h[..., None])[..., 0]
        dxx, dyy, dzz, dxy, dxz, dyz = beta[:, 1:].T
        T = np.stack([np.stack([dxx, dxy, dxz], 1), np.stack([dxy, dyy, dyz], 1), np.stack([dxz, dyz, dzz], 1)], 1)
        ev, evec = np.linalg.eigh(T)
        md = ev.mean(axis=1)
        fa = np.sqrt(1.5 * ((ev - md[:, None]) ** 2).sum(1) / np.maximum((ev ** 2).sum(1), 1e-12))
        return (fa, md, evec[:, :, 2])

    def noddi(self, y, v1):
        n = len(y)
        idx = np.empty(n, dtype=int)
        for s in range(0, n, 5000):
            idx[s:s + 5000] = np.argmax(np.abs(v1[s:s + 5000] @ DIRS.T), axis=1)
        out = np.full((n, 3), np.nan)
        for di in np.unique(idx):
            A, meta = self.dictionary(di)
            for r in np.where(idx == di)[0]:
                x, _ = nnls(A, y[r])
                xa, xi = (x[:-1], x[-1])
                s = xa.sum()
                if s > 0:
                    out[r] = (xa @ meta[:, 0] / s, xa @ meta[:, 1] / s, xi / (s + xi))
        return out

    def forward(self, icvf, odi, iso, mu):
        n = len(icvf)
        out = np.empty((n, len(self.bv)))
        for s in range(0, n, 1000):
            sl = slice(s, min(n, s + 1000))
            c2 = (mu[sl] @ SPHERE.T) ** 2
            kap = 1.0 / np.tan(odi[sl] * np.pi / 2)
            w = np.exp(kap[:, None] * (c2 - 1.0))
            w /= w.sum(1, keepdims=True)
            e_ic = w @ self.estick.T
            tau = (w * c2).sum(1)
            dperp = DPAR * (1 - icvf[sl])
            da, dr = (dperp + (DPAR - dperp) * tau, dperp + (DPAR - dperp) * (1 - tau) / 2)
            gm2 = (mu[sl] @ self.g.T) ** 2
            e_ec = np.exp(-self.bv[None, :] * (dr[:, None] + (da - dr)[:, None] * gm2))
            out[sl] = (1 - iso[sl])[:, None] * (icvf[sl][:, None] * e_ic + (1 - icvf[sl])[:, None] * e_ec) + iso[sl][:, None] * self.eiso[None, :]
        return out


def fit_methods(proto, S0, fats):
    res = {}
    S0 = S0.astype(np.float64)
    for name, fv in fats.items():
        S = S0 if fv is None else S0 - fv[:, None]
        b0n = S[:, proto.b0].mean(axis=1)
        ok = b0n > 0
        S = np.maximum(S, 0.001 * np.abs(b0n)[:, None])
        y = S / np.where(ok, b0n, 1.0)[:, None]
        fa, md, v1 = proto.dti(y)
        nd = proto.noddi(y, v1)
        met = {'FA': fa, 'MD': md, 'Vintra': nd[:, 0], 'ODI': nd[:, 1], 'Vcsf': nd[:, 2]}
        for key in met:
            met[key] = np.where(ok, met[key], np.nan)
        res[name] = met
    return res


SPHERE = fibonacci_sphere(6000)
_hemisphere = fibonacci_sphere(2000)
DIRS = _hemisphere[_hemisphere[:, 2] >= 0]
