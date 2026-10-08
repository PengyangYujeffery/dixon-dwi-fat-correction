import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
from dixon_dwi.core import calibrate, estimate_fat, shift_axis, subtract_fat
from dixon_dwi.cli import main
from dixon_dwi.models import Protocol, fit_methods, fibonacci_sphere
from dixon_dwi.comparators import fit_offset, physics_kernel

spec = importlib.util.spec_from_file_location("synthetic_demo", Path(__file__).parents[1] / "examples" / "synthetic_demo.py")
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)


class CoreTests(unittest.TestCase):
    def test_shift_direction_and_no_wrap(self):
        x = np.array([1., 0., 0., 0.])
        np.testing.assert_array_equal(shift_axis(x, 2), [0, 0, 1, 0])
        np.testing.assert_array_equal(shift_axis(x, -1), [0, 0, 0, 0])
        np.testing.assert_array_equal(shift_axis(x, 10), np.zeros(4))

    def test_known_signal_subtraction_including_b0(self):
        fat = np.arange(12.).reshape(4, 3)
        estimate = estimate_fat(fat, [-1, 1], [0.1, 0.2], scale=2)
        clean = np.full((4, 3, 5), 10.)
        np.testing.assert_allclose(subtract_fat(clean + estimate[..., None], estimate), clean)

    def test_artificial_recovery_on_separate_case(self):
        result = demo.run()
        self.assertAlmostEqual(result["kernel_mass"], 0.16, delta=0.01)
        self.assertAlmostEqual(result["centroid_pixels"], -4, delta=0.15)
        self.assertLess(result["corrected_rmse"], 0.1 * result["uncorrected_rmse"])

    def test_kernel_sign_follows_array_flip(self):
        original = [demo.make_case(seed)[0] for seed in (101, 102)]
        flipped = [{k: np.flip(v, 0) if np.ndim(v) else v for k, v in c.items()} for c in original]
        k = calibrate(original)
        f = calibrate(flipped)
        np.testing.assert_allclose(f["weights"], np.array(k["weights"])[::-1], atol=1e-8)

    def test_cli_roundtrip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            case, clean, observed = demo.make_case(201)
            np.savez(root / "input.npz", **case, dwi=observed)
            main(["calibrate", "--input", str(root / "input.npz"), "--pe-axis", "0", "--output", str(root / "kernel.json")])
            main(["correct", "--input", str(root / "input.npz"), "--kernel", str(root / "kernel.json"), "--output", str(root / "output.npz")])
            with np.load(root / "output.npz", allow_pickle=False) as data:
                self.assertEqual(data.files, ["corrected"])
                self.assertLess(np.sqrt(np.mean((data["corrected"] - clean) ** 2)), 0.1)
            with self.assertRaises(SystemExit):
                main(["calibrate", "--input", str(root / "input.npz"), "--pe-axis", "0", "--output", str(root / "kernel.json")])

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            estimate_fat(np.ones((3, 3)), [0], [-1])
        with self.assertRaises(ValueError):
            estimate_fat(np.ones((3, 3)), [0.5], [1])
        with self.assertRaises(ValueError):
            calibrate([])
        case = demo.make_case(101)[0]
        case["ff"][0, 0, 0] = np.nan
        with self.assertRaises(ValueError):
            calibrate([case])

    def test_physics_mass_and_polarity(self):
        shifts, negative = physics_kernel(20, 1.25, -1, 0.176)
        _, positive = physics_kernel(20, 1.25, 1, 0.176)
        self.assertAlmostEqual(negative.sum(), 0.176)
        self.assertLess(shifts @ negative, 0)
        np.testing.assert_allclose(negative, positive[::-1], atol=1e-14)

    def test_tensor_isotropic_reference(self):
        g = fibonacci_sphere(12)
        b = np.r_[0, np.full(12, 1000.), np.full(12, 2000.)]
        vectors = np.column_stack([np.zeros(3), g.T, g.T])
        proto = Protocol(b, vectors)
        fa, md, _ = proto.dti(np.exp(-b * 0.0015)[None, :])
        self.assertAlmostEqual(float(md[0]), 1.5, places=7)
        self.assertLess(float(fa[0]), 1e-7)


if __name__ == "__main__":
    unittest.main()
