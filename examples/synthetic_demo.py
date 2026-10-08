"""Run a small example built entirely from artificial signals."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter
from dixon_dwi.core import calibrate, estimate_fat, subtract_fat


def make_case(seed):
    rng = np.random.default_rng(seed)
    x, y, z = np.indices((96, 72, 4))
    fat_region = ((x >= 16) & (x < 30)) | ((x >= 66) & (x < 80))
    texture = gaussian_filter(rng.normal(size=x.shape), 1.0)
    fat = np.where(fat_region, 100.0 * (1 + 0.15 * texture), 0.0)
    water = np.where(fat_region, 0.0, 100.0 * (1 + 0.1 * texture))
    true_shifts, true_weights = [-5, -4, -3], [0.03, 0.10, 0.03]
    contamination = estimate_fat(fat, true_shifts, true_weights)
    bvalues = np.array([0, 500, 1000, 1500])
    clean = water[..., None] * np.exp(-bvalues * 0.0015)
    observed = clean + contamination[..., None]
    case = {"b0": observed[..., 0], "d1000": observed[..., 2], "fat": fat,
            "body": np.ones(x.shape, dtype=bool), "ff": fat_region.astype(float), "scale": 1.0}
    return case, clean, observed


def run():
    kernel = calibrate([make_case(101)[0], make_case(102)[0]])
    case, clean, observed = make_case(103)
    correction = estimate_fat(case["fat"], kernel["shifts"], kernel["weights"])
    corrected = subtract_fat(observed, correction)
    colocated = subtract_fat(observed, kernel["mass"] * case["fat"])
    rmse = lambda signal: float(np.sqrt(np.mean((signal - clean) ** 2)))
    result = {"example": "artificial signals only", "kernel_mass": kernel["mass"],
              "centroid_pixels": kernel["centroid_pixels"], "uncorrected_rmse": rmse(observed),
              "colocated_rmse": rmse(colocated), "corrected_rmse": rmse(corrected)}
    if not result["corrected_rmse"] < result["colocated_rmse"] < result["uncorrected_rmse"]:
        raise AssertionError("unexpected example result")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
