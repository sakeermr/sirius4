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
OUT = Path(sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-")
           else "examples/example_input.zip")


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


def make_deimos(n_feat: int = 6):
    """Synthetic tables in the DEIMoS paired-column layout."""
    ms1_rows, ms2_rows, feat = [], [], []
    for k in range(n_feat):
        pm = float(rng.uniform(150, 600)); rt = float(rng.uniform(1, 20))
        feat.append((k, pm, rt))
        for iso in range(3):
            ms1_rows.append(dict(controllerType=0, controllerNumber=1, scan=1000 + k,
                retention_time=round(rt, 4), mz=round(pm + iso * 1.00336, 5),
                intensity=round(1e6 * 0.6 ** iso, 1), persistence=round(1e5 * 0.6 ** iso, 1),
                mz_weighted=round(pm + iso * 1.00336, 5), retention_time_weighted=round(rt, 4)))
    for _ in range(400):  # background peaks that must not leak into a feature
        ms1_rows.append(dict(controllerType=0, controllerNumber=1, scan=9999,
            retention_time=round(float(rng.uniform(1, 20)), 4),
            mz=round(float(rng.uniform(150, 600)), 5), intensity=500.0,
            persistence=50.0, mz_weighted=0.0, retention_time_weighted=0.0))
    for k, pm, rt in feat:
        for j in range(int(rng.integers(6, 14))):
            ms2_rows.append(dict(index_ms1=k, mz_ms1=round(pm, 5),
                retention_time_ms1=round(rt, 4), intensity_ms1=1e6, persistence_ms1=1e5,
                index_ms2=j, mz_ms2=round(float(rng.uniform(50, pm)), 5),
                retention_time_ms2=round(rt, 4),
                intensity_ms2=round(float(rng.uniform(1e3, 5e5)), 1),
                persistence_ms2=1e4, retention_time_error=0.001))
    return pd.DataFrame(ms1_rows), pd.DataFrame(ms2_rows)


def main_deimos() -> None:
    tmp = OUT.parent / "_tmp_deimos"; tmp.mkdir(parents=True, exist_ok=True)
    stem = "20260401_SAMPLE_01"
    m1, m2 = make_deimos()
    m1.to_csv(tmp / f"{stem}_ms1_peaks.csv", index=False)
    m2.to_csv(tmp / f"{stem}_ms2_extracted.csv", index=False)
    m1.head(5).to_csv(tmp / f"{stem}_ms1_raw.csv", index=False)     # must be ignored
    m2.head(5).to_csv(tmp / f"{stem}_ms2_peaks.csv", index=False)   # must be ignored
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(tmp.iterdir()):
            zf.write(p, arcname=p.name)
    for p in tmp.iterdir():
        p.unlink()
    tmp.rmdir()
    print(f"wrote {OUT} (DEIMoS layout)")


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
    if "--deimos" in sys.argv:
        main_deimos()
    else:
        main()
