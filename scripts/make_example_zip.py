#!/usr/bin/env python3
"""
Creates a small synthetic input archive with the same shape as the real data:

    example_input.zip
      ms1_peaks_batch1.csv
      ms1_peaks_batch2.h5
      ms2_extracted_batch1.csv
      ms2_extracted_batch2.h5

Used by the smoke-test workflow so the pipeline can be validated without
uploading real data. Run:  python scripts/make_example_zip.py examples/example_input.zip
"""
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "examples/example_input.zip")


def make(prefix: str, n_feat: int, start: int):
    ms1_rows, ms2_rows = [], []
    for k in range(n_feat):
        fid = f"FT{start + k:04d}"
        pmz = float(rng.uniform(150, 600))
        rt = float(rng.uniform(30, 900))
        # MS1 isotope pattern
        for iso in range(3):
            ms1_rows.append({
                "feature_id": fid,
                "mz": round(pmz + iso * 1.00336, 5),
                "intensity": round(1e6 * (0.6 ** iso), 2),
                "rt": round(rt, 2),
                "adduct": "[M+H]+",
            })
        # MS2 fragments at two collision energies
        for ce in (20, 40):
            for _ in range(rng.integers(5, 12)):
                ms2_rows.append({
                    "feature_id": fid,
                    "precursor_mz": round(pmz, 5),
                    "mz": round(float(rng.uniform(50, pmz)), 5),
                    "intensity": round(float(rng.uniform(1e3, 5e5)), 2),
                    "rt": round(rt, 2),
                    "collision_energy": ce,
                    "charge": 1,
                })
    return pd.DataFrame(ms1_rows), pd.DataFrame(ms2_rows)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.parent / "_tmp_example"
    tmp.mkdir(exist_ok=True)

    ms1a, ms2a = make("b1", 4, 1)
    ms1b, ms2b = make("b2", 3, 101)

    ms1a.to_csv(tmp / "ms1_peaks_batch1.csv", index=False)
    ms2a.to_csv(tmp / "ms2_extracted_batch1.csv", index=False)
    ms1b.to_hdf(tmp / "ms1_peaks_batch2.h5", key="ms1_peaks", mode="w")
    ms2b.to_hdf(tmp / "ms2_extracted_batch2.h5", key="ms2_extracted", mode="w")

    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(tmp.iterdir()):
            zf.write(p, arcname=p.name)
    for p in tmp.iterdir():
        p.unlink()
    tmp.rmdir()
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
