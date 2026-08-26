#!/usr/bin/env python3
"""
Validate generated .ms files. Used by the smoke-test workflow.

    python scripts/check_ms_files.py <dir> --min 6 --require-isotopes

Checks that the >ms1peaks block really is an isotope envelope (peaks spaced by
~1.0034 Da) and not a cloud of co-eluting noise, which is what a naive m/z
window produces on a large peak list.
"""
import argparse
import glob
import sys
from pathlib import Path

NEUTRON = 1.0033548


def ms1_peaks(path: Path) -> list[float]:
    mzs, inside = [], False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">ms1peaks"):
            inside = True
            continue
        if line.startswith(">"):
            inside = False
            continue
        if inside and line.strip():
            try:
                mzs.append(float(line.split()[0]))
            except ValueError:
                pass
    return mzs


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("directory")
    ap.add_argument("--min", type=int, default=1, help="minimum number of .ms files")
    ap.add_argument("--require-isotopes", action="store_true",
                    help="at least one file must have a 2+ peak isotope envelope")
    ap.add_argument("--tol", type=float, default=0.02, help="allowed deviation from 1.0034 Da")
    args = ap.parse_args()

    files = sorted(glob.glob(f"{args.directory}/*.ms"))
    print(f"found {len(files)} .ms files in {args.directory}")
    if len(files) < args.min:
        print(f"FAIL: expected at least {args.min}")
        return 1

    bad, with_env = 0, 0
    for f in files:
        mzs = ms1_peaks(Path(f))
        if len(mzs) > 1:
            with_env += 1
            gaps = [round(b - a, 4) for a, b in zip(mzs, mzs[1:])]
            off = [g for g in gaps if abs(g - NEUTRON) > args.tol]
            if off:
                print(f"FAIL: {Path(f).name} MS1 peaks are not an isotope pattern; gaps={gaps}")
                bad += 1

    print(f"{with_env}/{len(files)} files have a multi-peak MS1 envelope, {bad} malformed")
    if bad:
        return 1
    if args.require_isotopes and with_env == 0:
        print("FAIL: no file got a real isotope envelope")
        return 1

    sample = Path(files[0])
    print(f"--- sample: {sample.name} ---")
    print("\n".join(sample.read_text(encoding="utf-8").splitlines()[:16]))
    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
