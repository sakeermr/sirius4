#!/usr/bin/env python3
"""
convert_to_ms.py
================

Takes a ZIP archive (or an already-extracted folder) that contains

    *ms1_peaks*.csv / .h5
    *ms2_extracted*.csv / .h5

and converts every feature it finds into a SIRIUS `.ms` file.

SIRIUS cannot read CSV or HDF5, so this step is mandatory before running
the `sirius` CLI inside the rformassspectrometry/rusirius container.

Output layout
-------------
    <out>/ms_files/FT0001.ms
    <out>/ms_files/FT0002.ms
    ...
    <out>/manifest.csv          one row per written compound
    <out>/conversion_report.txt human readable log

Usage
-----
    python scripts/convert_to_ms.py --zip data/input.zip --out work
    python scripts/convert_to_ms.py --dir some/folder   --out work \
        --adduct "[M+H]+" --polarity positive
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# column-name aliases -> canonical name
# add your own spellings here (or in config/columns.yaml) if a file is missed
# --------------------------------------------------------------------------
ALIASES: dict[str, list[str]] = {
    "feature_id": [
        "feature_id", "featureid", "feature", "id", "compound_id", "compoundid",
        "group_id", "groupid", "cluster_id", "clusterid", "scan_id", "scanid",
        "spectrum_id", "spectrumid", "index", "idx", "name", "title", "key",
    ],
    "mz": ["mz", "m/z", "m_z", "mass", "mzs", "peak_mz", "mz_value", "fragment_mz", "mass_to_charge"],
    "intensity": [
        "intensity", "intensities", "int", "i", "abundance", "abundances",
        "height", "area", "peak_intensity", "counts", "signal",
    ],
    "rt": ["rt", "retention_time", "retentiontime", "rtime", "rt_sec", "rt_seconds", "rt_min", "rtmed", "rt_apex"],
    "precursor_mz": [
        "precursor_mz", "precursormz", "precursor", "parent_mass", "parentmass",
        "parent_mz", "premz", "pepmass", "mz_precursor", "precursor_mass",
        "isolation_mz", "target_mz",
    ],
    "charge": ["charge", "z", "precursor_charge"],
    "adduct": ["adduct", "ionization", "ion", "ion_type", "adduct_type", "precursor_type", "ionmode_adduct"],
    "collision_energy": [
        "collision_energy", "collisionenergy", "ce", "energy", "hcd", "hcd_energy",
        "nce", "normalized_collision_energy", "collision",
    ],
    "formula": ["formula", "molecular_formula", "molecularformula", "mf", "sum_formula"],
    "compound_name": ["compound_name", "compound", "molecule", "annotation", "label"],
    "polarity": ["polarity", "ion_mode", "ionmode", "mode", "msmode"],
}

MS1_PATTERNS = ["*ms1_peak*", "*ms1peak*", "*ms1*"]
MS2_PATTERNS = ["*ms2_extract*", "*ms2extract*", "*ms2*"]
TABLE_EXT = {".csv", ".tsv", ".txt", ".h5", ".hdf5", ".hdf", ".he5"}

LOG: list[str] = []


def log(msg: str) -> None:
    print(msg, flush=True)
    LOG.append(str(msg))


# --------------------------------------------------------------------------
# I/O helpers
# --------------------------------------------------------------------------
def unzip(zip_path: Path, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        # protect against zip-slip
        for member in zf.namelist():
            target = (dest / member).resolve()
            if not str(target).startswith(str(dest.resolve())):
                raise RuntimeError(f"Unsafe path in archive: {member}")
        zf.extractall(dest)
    log(f"[unzip] extracted {zip_path.name} -> {dest}")
    return dest


def find_tables(root: Path, patterns: list[str]) -> list[Path]:
    """Find data files anywhere under root matching any of the glob patterns."""
    hits: list[Path] = []
    for pat in patterns:
        for p in root.rglob("*"):
            if not p.is_file():
                continue
            if p.suffix.lower() not in TABLE_EXT:
                continue
            if p.name.startswith("._") or "__MACOSX" in str(p):
                continue
            if Path(p.name.lower()).match(pat):
                hits.append(p)
        if hits:  # first pattern that matches wins (most specific first)
            break
    return sorted(set(hits))


def read_csv_like(path: Path) -> list[pd.DataFrame]:
    sep = "\t" if path.suffix.lower() in {".tsv"} else None
    df = pd.read_csv(path, sep=sep, engine="python")
    return [df]


def read_h5(path: Path) -> list[pd.DataFrame]:
    """Read an HDF5 file. Tries pandas/PyTables first, falls back to h5py."""
    frames: list[pd.DataFrame] = []
    try:
        with pd.HDFStore(str(path), mode="r") as store:
            for key in store.keys():
                obj = store.get(key)
                if isinstance(obj, pd.Series):
                    obj = obj.to_frame()
                if isinstance(obj, pd.DataFrame) and len(obj):
                    obj = obj.reset_index()
                    obj.attrs["h5_key"] = key
                    frames.append(obj)
        if frames:
            log(f"[read] {path.name}: {len(frames)} pandas/HDF5 table(s)")
            return frames
    except Exception as exc:  # not a pandas HDFStore -> generic HDF5
        log(f"[read] {path.name}: not a pandas HDFStore ({type(exc).__name__}), using h5py")

    import h5py

    columns: dict[str, np.ndarray] = {}

    def visit(name: str, obj) -> None:
        if not isinstance(obj, h5py.Dataset):
            return
        data = obj[()]
        short = name.split("/")[-1]
        if isinstance(data, np.ndarray) and data.dtype.names:  # compound/structured
            frames.append(pd.DataFrame({n: _flat(data[n]) for n in data.dtype.names}))
        elif isinstance(data, np.ndarray) and data.ndim == 2 and data.shape[1] in (2, 3):
            names = ["mz", "intensity", "rt"][: data.shape[1]]
            frames.append(pd.DataFrame(data, columns=names))
        elif isinstance(data, np.ndarray) and data.ndim == 1:
            columns[short] = _flat(data)

    with h5py.File(path, "r") as f:
        f.visititems(visit)

    if columns:
        n = max(len(v) for v in columns.values())
        aligned = {k: v for k, v in columns.items() if len(v) == n}
        if aligned:
            frames.append(pd.DataFrame(aligned))
    log(f"[read] {path.name}: {len(frames)} table(s) via h5py")
    return frames


def _flat(arr: np.ndarray) -> np.ndarray:
    if arr.dtype.kind in {"S", "O"}:
        return np.array([x.decode() if isinstance(x, bytes) else x for x in arr])
    return arr


def load_any(path: Path) -> list[pd.DataFrame]:
    if path.suffix.lower() in {".h5", ".hdf5", ".hdf", ".he5"}:
        return read_h5(path)
    return read_csv_like(path)


# --------------------------------------------------------------------------
# normalisation
# --------------------------------------------------------------------------
def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def canonicalise(df: pd.DataFrame, extra_aliases: dict | None = None) -> pd.DataFrame:
    aliases = {k: list(v) for k, v in ALIASES.items()}
    for k, v in (extra_aliases or {}).items():
        aliases.setdefault(k, [])
        aliases[k] = list(v) + aliases[k]

    lookup = {}
    for canon, names in aliases.items():
        for n in names:
            lookup.setdefault(_norm(n), canon)

    rename: dict[str, str] = {}
    used: set[str] = set()
    for col in df.columns:
        canon = lookup.get(_norm(col))
        if canon and canon not in used:
            rename[col] = canon
            used.add(canon)
    out = df.rename(columns=rename).copy()
    out.columns = [str(c) for c in out.columns]
    return out


def to_float(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


# --------------------------------------------------------------------------
# building compounds
# --------------------------------------------------------------------------
class Compound:
    def __init__(self, fid: str):
        self.fid = str(fid)
        self.ms1: list[tuple[float, float]] = []
        self.ms2: dict[str, list[tuple[float, float]]] = defaultdict(list)
        self.precursor_mz: float | None = None
        self.rt: float | None = None
        self.charge: int | None = None
        self.adduct: str | None = None
        self.formula: str | None = None
        self.name: str | None = None

    def parent_mass(self) -> float | None:
        if self.precursor_mz and self.precursor_mz > 0:
            return float(self.precursor_mz)
        if self.ms1:
            return float(max(self.ms1, key=lambda p: p[1])[0])
        return None

    def n_ms2(self) -> int:
        return sum(len(v) for v in self.ms2.values())


def safe_id(x) -> str:
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", str(x)).strip("_")
    return s or "unknown"


def collect_ms1(frames: list[pd.DataFrame], comps: dict[str, Compound], min_int: float) -> None:
    for df in frames:
        df = canonicalise(df)
        if "mz" not in df.columns:
            log(f"[warn] MS1 table without an m/z column, columns={list(df.columns)[:12]} -- skipped")
            continue
        df["mz"] = to_float(df["mz"])
        df["intensity"] = to_float(df["intensity"]) if "intensity" in df.columns else 1.0
        df = df.dropna(subset=["mz"])
        df = df[df["intensity"].fillna(0) >= min_int]
        if "feature_id" not in df.columns:
            if "precursor_mz" in df.columns and to_float(df["precursor_mz"]).notna().any():
                df["feature_id"] = to_float(df["precursor_mz"]).round(4).astype(str)
                log("[info] MS1 table had no feature id -> grouping by rounded precursor m/z "
                    "(keeps MS1 and MS2 of the same feature together)")
            else:
                df["feature_id"] = [f"F{i+1:05d}" for i in range(len(df))]
                log("[info] MS1 table had no feature id column -> one feature per row")
        for fid, grp in df.groupby(df["feature_id"].map(safe_id)):
            c = comps.setdefault(fid, Compound(fid))
            c.ms1.extend(zip(grp["mz"].astype(float), grp["intensity"].fillna(1.0).astype(float)))
            if "rt" in grp.columns and c.rt is None:
                v = to_float(grp["rt"]).dropna()
                if len(v):
                    c.rt = float(v.iloc[0])
            if "precursor_mz" in grp.columns and c.precursor_mz is None:
                v = to_float(grp["precursor_mz"]).dropna()
                if len(v):
                    c.precursor_mz = float(v.iloc[0])
            for attr, col in (("adduct", "adduct"), ("formula", "formula"), ("name", "compound_name")):
                if col in grp.columns and getattr(c, attr) is None:
                    v = grp[col].dropna()
                    if len(v):
                        setattr(c, attr, str(v.iloc[0]))


def collect_ms2(frames: list[pd.DataFrame], comps: dict[str, Compound], min_int: float) -> None:
    for df in frames:
        df = canonicalise(df)
        if "mz" not in df.columns:
            log(f"[warn] MS2 table without an m/z column, columns={list(df.columns)[:12]} -- skipped")
            continue
        df["mz"] = to_float(df["mz"])
        df["intensity"] = to_float(df["intensity"]) if "intensity" in df.columns else 1.0
        df = df.dropna(subset=["mz"])
        df = df[df["intensity"].fillna(0) >= min_int]
        if "feature_id" not in df.columns:
            if "precursor_mz" in df.columns:
                df["feature_id"] = to_float(df["precursor_mz"]).round(4).astype(str)
                log("[info] MS2 table had no feature id -> grouping by rounded precursor m/z")
            else:
                df["feature_id"] = "F00001"
        for fid, grp in df.groupby(df["feature_id"].map(safe_id)):
            c = comps.setdefault(fid, Compound(fid))
            if "precursor_mz" in grp.columns:
                v = to_float(grp["precursor_mz"]).dropna()
                if len(v):
                    c.precursor_mz = float(v.iloc[0])
            if "rt" in grp.columns and c.rt is None:
                v = to_float(grp["rt"]).dropna()
                if len(v):
                    c.rt = float(v.iloc[0])
            if "charge" in grp.columns and c.charge is None:
                v = to_float(grp["charge"]).dropna()
                if len(v):
                    c.charge = int(v.iloc[0])
            for attr, col in (("adduct", "adduct"), ("formula", "formula"), ("name", "compound_name")):
                if col in grp.columns and getattr(c, attr) is None:
                    v = grp[col].dropna()
                    if len(v):
                        setattr(c, attr, str(v.iloc[0]))
            if "collision_energy" in grp.columns:
                for ce, sub in grp.groupby(grp["collision_energy"].fillna("NA").astype(str)):
                    c.ms2[str(ce)].extend(
                        zip(sub["mz"].astype(float), sub["intensity"].fillna(1.0).astype(float))
                    )
            else:
                c.ms2["NA"].extend(
                    zip(grp["mz"].astype(float), grp["intensity"].fillna(1.0).astype(float))
                )


# --------------------------------------------------------------------------
# .ms writer
# --------------------------------------------------------------------------
def write_ms(c: Compound, out_dir: Path, args) -> Path | None:
    pm = c.parent_mass()
    if pm is None:
        return None
    adduct = c.adduct or args.adduct
    charge = c.charge if c.charge else (1 if args.polarity == "positive" else -1)
    if adduct and adduct.strip().endswith("-"):
        charge = -abs(charge)

    lines: list[str] = []
    lines.append(f">compound {c.name or c.fid}")
    lines.append(f">parentmass {pm:.5f}")
    if adduct:
        lines.append(f">ionization {adduct}")
    lines.append(f">charge {charge}")
    if c.formula:
        lines.append(f">formula {c.formula}")
    if c.rt is not None and np.isfinite(c.rt):
        lines.append(f">rt {c.rt:.3f}")
    if args.instrument:
        lines.append(f">instrumentation {args.instrument}")
    lines.append("")

    if c.ms1:
        lines.append(">ms1peaks")
        for mz, it in sorted(dedupe(c.ms1)):
            lines.append(f"{mz:.5f} {it:.4f}")
        lines.append("")

    for ce, peaks in c.ms2.items():
        if not peaks:
            continue
        header = ">ms2peaks" if ce in {"NA", "nan", ""} else f">collision {strip_ce(ce)}"
        lines.append(header)
        for mz, it in sorted(dedupe(peaks)):
            lines.append(f"{mz:.5f} {it:.4f}")
        lines.append("")

    path = out_dir / f"{safe_id(c.fid)}.ms"
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return path


def strip_ce(ce: str) -> str:
    m = re.search(r"[-+]?\d*\.?\d+", str(ce))
    return m.group(0) if m else str(ce)


def dedupe(peaks: list[tuple[float, float]]) -> list[tuple[float, float]]:
    agg: dict[float, float] = {}
    for mz, it in peaks:
        if not np.isfinite(mz) or mz <= 0:
            continue
        key = round(float(mz), 5)
        agg[key] = max(agg.get(key, 0.0), float(it) if np.isfinite(it) else 0.0)
    return sorted(agg.items())


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description="Convert ms1_peaks / ms2_extracted tables to SIRIUS .ms files")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--zip", type=Path, help="input ZIP archive")
    src.add_argument("--dir", type=Path, help="already extracted folder")
    ap.add_argument("--out", type=Path, default=Path("work"), help="output folder")
    ap.add_argument("--adduct", default="[M+H]+", help="default adduct if the tables do not carry one")
    ap.add_argument("--polarity", choices=["positive", "negative"], default="positive")
    ap.add_argument("--instrument", default="", help='e.g. "Orbitrap" or "Q-TOF"')
    ap.add_argument("--min-intensity", type=float, default=0.0)
    ap.add_argument("--min-ms2-peaks", type=int, default=1, help="skip compounds with fewer MS2 peaks")
    ap.add_argument("--require-ms2", action="store_true", help="only write compounds that have MS2")
    ap.add_argument("--max-compounds", type=int, default=0, help="0 = no limit (useful for smoke tests)")
    args = ap.parse_args()

    out = args.out
    ms_dir = out / "ms_files"
    ms_dir.mkdir(parents=True, exist_ok=True)

    root = unzip(args.zip, out / "extracted") if args.zip else args.dir
    root = Path(root)
    if not root.exists():
        log(f"[error] input not found: {root}")
        return 2

    # if the zip contained a single top-level folder, dive into it
    kids = [p for p in root.iterdir() if not p.name.startswith("__")]
    if len(kids) == 1 and kids[0].is_dir():
        log(f"[info] descending into single top-level folder: {kids[0].name}")

    ms1_files = find_tables(root, MS1_PATTERNS)
    ms2_files = find_tables(root, MS2_PATTERNS)
    ms2_files = [p for p in ms2_files if p not in set(ms1_files)]
    ms1_files = [p for p in ms1_files if "ms2" not in p.name.lower()]

    log(f"[scan] MS1 files ({len(ms1_files)}): {[p.name for p in ms1_files]}")
    log(f"[scan] MS2 files ({len(ms2_files)}): {[p.name for p in ms2_files]}")
    if not ms1_files and not ms2_files:
        log("[error] no ms1_peaks / ms2_extracted tables found in the archive")
        return 3

    comps: dict[str, Compound] = {}
    for p in ms1_files:
        try:
            collect_ms1(load_any(p), comps, args.min_intensity)
        except Exception as exc:
            log(f"[warn] failed to read {p}: {type(exc).__name__}: {exc}")
    for p in ms2_files:
        try:
            collect_ms2(load_any(p), comps, args.min_intensity)
        except Exception as exc:
            log(f"[warn] failed to read {p}: {type(exc).__name__}: {exc}")

    log(f"[build] {len(comps)} candidate features")

    written, skipped = [], 0
    for fid, c in sorted(comps.items()):
        if args.require_ms2 and c.n_ms2() < max(1, args.min_ms2_peaks):
            skipped += 1
            continue
        p = write_ms(c, ms_dir, args)
        if p is None:
            skipped += 1
            continue
        written.append((fid, p, c))
        if args.max_compounds and len(written) >= args.max_compounds:
            break

    with (out / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["feature_id", "ms_file", "parent_mass", "rt", "adduct", "n_ms1_peaks", "n_ms2_peaks"])
        for fid, p, c in written:
            w.writerow([
                fid, p.name, f"{c.parent_mass():.5f}",
                "" if c.rt is None else f"{c.rt:.3f}",
                c.adduct or args.adduct, len(dedupe(c.ms1)), c.n_ms2(),
            ])

    log(f"[done] wrote {len(written)} .ms files to {ms_dir} (skipped {skipped})")
    (out / "conversion_report.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    (out / "conversion_summary.json").write_text(
        json.dumps({"n_written": len(written), "n_skipped": skipped,
                    "ms1_files": [p.name for p in ms1_files],
                    "ms2_files": [p.name for p in ms2_files]}, indent=2),
        encoding="utf-8",
    )
    return 0 if written else 4


if __name__ == "__main__":
    sys.exit(main())
