"""Check the tracked release against an explicit file list."""
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
ALLOWED = set((ROOT / "release_files.txt").read_text(encoding="utf-8").splitlines())
PATTERNS = {
    "absolute local path": r"(?i)(?:[A-Z]:[\\/](?:Users|home)[\\/]|/home/[^\s/]+/)",
    "case code": r"\b(?:HC|KP)\d{2,}\b",
    "participant field": r"(?i)(?:PatientName|PatientID|PatientBirthDate|AccessionNumber|a_by_subject|_included_sids|_skipped_sids|per_subject_A_PE)",
    "email": r"[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}",
    "private key": r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----",
    "credential": r"(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|AKIA[A-Z0-9]{16})",
    "embedded binary": r"data:[\w/+.-]+;base64,",
}


def files_to_check():
    if (ROOT / ".git").exists():
        raw = subprocess.check_output(["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"], cwd=ROOT)
        return {p for p in raw.decode("utf-8").split("\0") if p}
    files = set()
    for path in ROOT.rglob("*"):
        relative = path.relative_to(ROOT)
        if any(part in {".git", ".venv", "__pycache__", "build", "dist"} or part.endswith(".egg-info") for part in relative.parts):
            continue
        if path.is_file():
            files.add(relative.as_posix())
    return files


def main():
    files = files_to_check()
    errors = []
    for name in sorted(files - ALLOWED):
        errors.append(f"not allowlisted: {name}")
    for name in sorted(ALLOWED - files):
        errors.append(f"missing: {name}")
    for name in sorted(files & ALLOWED):
        path = ROOT / name
        if path.is_symlink() or not path.resolve().is_relative_to(ROOT):
            errors.append(f"unsafe path: {name}")
            continue
        if path.stat().st_size > 250_000:
            errors.append(f"oversized file: {name}")
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeError, OSError):
            errors.append(f"non-text file: {name}")
            continue
        if name == "tools/check_release.py":
            continue
        for label, pattern in PATTERNS.items():
            if re.search(pattern, text):
                errors.append(f"{label}: {name}")
    if errors:
        print("Release check failed:")
        print("\n".join(errors))
        return 1
    print(f"Release check passed: {len(files)} allowlisted text files.")
    print("No dataset or binary asset is permitted. Review contents before publication.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
