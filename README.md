# sirius4

Automated SIRIUS annotation pipeline on GitHub Actions.

You upload a **ZIP** containing `ms1_peaks*` and `ms2_extracted*` tables (CSV **and/or** HDF5).
The pipeline unzips it, converts every feature into a SIRIUS **`.ms`** file, and runs the
SIRIUS CLI inside the Docker image `rformassspectrometry/rusirius`.
Results come back as downloadable workflow artifacts.

```
input.zip ──► convert_to_ms.py ──► work/ms_files/*.ms ──► docker run rusirius ──► work/project/
   csv/h5           (python)          SIRIUS format          (sirius CLI)         + summaries
```

---

## 1. Create the repo

```bash
# with the GitHub CLI
gh repo create sirius4 --public --source=. --remote=origin --push

# or manually
git init && git add . && git commit -m "sirius4 pipeline"
git branch -M main
git remote add origin https://github.com/<YOUR-USER>/sirius4.git
git push -u origin main
```

## 2. Add your SIRIUS account (required for CSI:FingerID / CANOPUS)

SIRIUS 5.8+ needs a free account for the web-service steps (`fingerprint`, `structure`, `canopus`).
The `formula` step works without one.

**Settings → Secrets and variables → Actions → New repository secret**

| Secret | Value |
| --- | --- |
| `SIRIUS_USER` | your SIRIUS/Bright Giant e-mail |
| `SIRIUS_PASSWORD` | your password |

## 3. Run it

**Option A — commit the zip**

```bash
cp /path/to/my_data.zip data/input.zip
git add data/input.zip && git commit -m "add data" && git push
```
The workflow fires automatically. (GitHub's hard file limit is 100 MB; for bigger data see
*Large inputs* below.)

**Option B — Actions tab → “SIRIUS pipeline” → Run workflow**, and set the inputs:

| Input | Default | Meaning |
| --- | --- | --- |
| `input_zip` | `data/input.zip` | path to the archive in the repo |
| `adduct` | `[M+H]+` | fallback adduct if your tables don't have one |
| `polarity` | `positive` | sets the charge sign |
| `instrument` | `orbitrap` | `orbitrap` or `qtof` |
| `ppm_max` | `10` | MS1/MS2 mass accuracy |
| `tools` | `formula fingerprint structure canopus` | SIRIUS subcommands |
| `max_compounds` | `0` | limit for a quick test run |

## 4. Get the results

Open the finished run → **Artifacts**:

* `ms-files` — the generated `.ms` files, `manifest.csv`, `conversion_report.txt`
* `sirius-results` — the full SIRIUS project directory plus `summaries/*.tsv`

---

## Repository layout

```
sirius4/
├── .github/workflows/
│   ├── sirius.yml          full pipeline: convert -> SIRIUS
│   └── smoke-test.yml      converter test on synthetic data (no Docker)
├── scripts/
│   ├── convert_to_ms.py    CSV/H5 -> SIRIUS .ms converter
│   ├── run_sirius.sh       runs the sirius CLI in the rusirius container
│   └── make_example_zip.py builds a fake input zip for testing
├── config/columns.yaml     optional extra column-name spellings
├── data/                   put your input .zip here
├── examples/               sample .ms output
├── requirements.txt
└── Makefile
```

## Run locally

```bash
make deps
make convert IN=data/input.zip      # -> work/ms_files/*.ms
make sirius                          # needs Docker
```

---

## DEIMoS input

DEIMoS (`pnnl/deimos`) writes `ms2_extracted` tables where each row is one
**precursor/fragment pair**, not a plain peak list:

| index_ms1 | mz_ms1 | retention_time_ms1 | index_ms2 | mz_ms2 | intensity_ms2 | ... |
|---|---|---|---|---|---|---|

`index_ms1` is the feature id, `mz_ms1` the precursor, `mz_ms2` the fragment.
The `ms1_peaks` table has **no** feature id - it is a flat peak list - so the
isotope pattern of each feature is recovered by searching that list in a window
around the precursor (`-0.5` to `+4.5` Da, and +/- `--iso-rt-tol` in retention time).

DEIMoS writes `ms2_extracted` in one of two shapes, both handled automatically:

* **one row per feature**, with the whole fragment spectrum inside a single
  cell as a list — `mz_ms2 = "[283.1687, 223.1113, ...]"`. This is what
  `deimos.alignment`/`deimos.peakpick` normally produce.
* **one row per precursor/fragment pair**, with scalar `mz_ms2` values.

Two knobs matter:

* `--iso-ppm` / `--n-isotopes` - the MS1 envelope is built by matching each
  expected isotope position (M, M+1.00335, M+2.0067 ...) within this ppm
  tolerance and taking the most intense peak. A plain m/z window does not work
  on a million-peak list: it scoops up every co-eluting ion near the precursor
  and hands SIRIUS something that is not an isotope pattern at all.
* `--iso-rt-tol` - **in the same unit as your `retention_time` column.** DEIMoS
  usually writes minutes, so `0.05` is 3 s. If your column is in seconds, use
  something like `3`. Too wide pulls in co-eluting noise; too narrow leaves the
  MS1 block empty. Check `conversion_report.txt`: it reports how many features
  got isotope peaks attached.
* `--prefer-format` - DEIMoS exports every table as **both** `.csv` and `.h5`.
  The converter keeps only one (default `h5`) so the data isn't loaded twice.

Files named `*_ms1_raw*`, `*_ms2_raw*` and `*_ms2_peaks*` are ignored - only
`ms1_peaks` and `ms2_extracted` are read.

## Inputs larger than 100 MB

GitHub rejects any file over 100 MB, so a big archive cannot be committed.
Upload it as a **release asset** instead - the workflow falls back to that
automatically when `data/` has no zip:

```bash
gh release create data --title "input data" --notes "SIRIUS input"
gh release upload data "path/to/deimos-results.zip"
```

Then run the workflow normally; it downloads the asset at the start of the job.
Assets can be up to 2 GB.

## Input format

The converter is deliberately forgiving. It finds any file whose name contains
`ms1_peak` / `ms2_extract` (any depth inside the zip, `.csv`, `.tsv`, `.h5`, `.hdf5`)
and maps columns by name, case- and separator-insensitive:

| Meaning | Recognised column names |
| --- | --- |
| feature id | `feature_id`, `id`, `compound_id`, `group_id`, `scan_id`, `name`, … |
| m/z | `mz`, `m/z`, `mass`, `peak_mz`, `fragment_mz`, … |
| intensity | `intensity`, `int`, `abundance`, `height`, `area`, `counts`, … |
| retention time | `rt`, `retention_time`, `rtmed`, `rt_sec`, … |
| precursor | `precursor_mz`, `parentmass`, `pepmass`, `premz`, … |
| collision energy | `collision_energy`, `ce`, `nce`, `hcd`, … |
| adduct | `adduct`, `ionization`, `precursor_type`, … |

If one of your columns isn't picked up, check `work/conversion_report.txt` and add the
spelling to `ALIASES` in `scripts/convert_to_ms.py` (or `config/columns.yaml`).

**HDF5**: read via pandas/PyTables first; if the file isn't a pandas `HDFStore`, it falls back
to `h5py` and handles structured arrays, `N×2` (m/z, intensity) matrices, and parallel 1-D
datasets of equal length.

**Rows are grouped per feature.** If there is no feature-id column, features are grouped by
rounded precursor m/z; if there's no precursor either, each MS1 row becomes its own feature.

## Output format (`.ms`)

```
>compound FT0001
>parentmass 380.31973
>ionization [M+H]+
>charge 1
>rt 856.900

>ms1peaks
380.31973 1000000.0000
381.32309 600000.0000

>collision 20
97.61876 474376.0700
...

>collision 40
...
```

Spectra without a collision energy are written under `>ms2peaks`.

---

## Notes and gotchas

* **SIRIUS 5 vs 6 subcommand names.** The tool chain (`formula fingerprint
  structure canopus`) uses SIRIUS 5 names. SIRIUS 6 renamed several of them, so
  `run_sirius.sh` reads the image's own `--help`, remaps each requested tool to
  the name that build actually knows, and skips any that do not exist - logging
  every substitution. The resolved chain is printed as
  `==> resolved tool chain: ...` before the run.
* **The sirius binary path.** `run_sirius.sh` auto-discovers the executable inside the image
  (`command -v sirius`, then a `find` under `/opt`, `/usr/local`, …). If the image lays it out
  unusually, set `SIRIUS_BIN` explicitly:
  ```bash
  docker run --rm -it --entrypoint /bin/bash rformassspectrometry/rusirius
  SIRIUS_BIN=/opt/sirius/bin/sirius ./scripts/run_sirius.sh work/ms_files work/project
  ```
* **Runtime.** GitHub-hosted runners cap a job at 6 h and give 2 CPU / 7 GB RAM. `structure`
  and `canopus` on thousands of compounds will exceed that — use `max_compounds`, split the
  input, or move to a self-hosted runner (`runs-on: self-hosted`).
* **Large inputs (>100 MB).** Don't commit the zip. Either use
  `workflow_dispatch` + a download step (S3/Zenodo/release asset), or Git LFS.
* **Validate first.** The *Smoke test* workflow runs the converter against synthetic data on
  every push, so a broken parser is caught before you burn an hour of SIRIUS compute.

## License

MIT
