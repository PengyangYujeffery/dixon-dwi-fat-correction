"""Local array-file interface."""
import argparse
import json
from pathlib import Path
import numpy as np
from .core import calibrate, estimate_fat, subtract_fat


def load_case(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: np.array(data[key]) for key in ("b0", "d1000", "fat", "body", "ff", "scale")}


def main(argv=None):
    parser = argparse.ArgumentParser(description="Dixon-guided residual-fat correction")
    commands = parser.add_subparsers(dest="command", required=True)
    cal = commands.add_parser("calibrate", help="fit a kernel on prepared calibration arrays")
    cal.add_argument("--input", type=Path, action="append", required=True)
    cal.add_argument("--output", type=Path, required=True)
    cal.add_argument("--pe-axis", type=int, choices=[0, 1], required=True)
    cal.add_argument("--slices", type=int, nargs="+")
    correct = commands.add_parser("correct", help="correct prepared DWI locally")
    correct.add_argument("--input", type=Path, required=True)
    correct.add_argument("--kernel", type=Path, required=True)
    correct.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("output already exists; choose a new file")
    if args.command == "calibrate":
        if args.output.suffix.lower() != ".json":
            parser.error("calibration output must end in .json")
        kernel = calibrate((load_case(p) for p in args.input), pe_axis=args.pe_axis, slices=args.slices)
        with args.output.open("x", encoding="utf-8") as handle:
            json.dump(kernel, handle, indent=2, allow_nan=False)
            handle.write("\n")
    else:
        if args.output.suffix.lower() != ".npz":
            parser.error("correction output must end in .npz")
        kernel = json.loads(args.kernel.read_text(encoding="utf-8"))
        with np.load(args.input, allow_pickle=False) as data:
            fat = estimate_fat(data["fat"], kernel["shifts"], kernel["weights"],
                               float(data["scale"]), kernel["pe_axis"])
            corrected = subtract_fat(data["dwi"], fat)
        with args.output.open("xb") as handle:
            np.savez_compressed(handle, corrected=corrected)
    print("Done. Outputs remain local; review them before sharing.")


if __name__ == "__main__":
    main()
