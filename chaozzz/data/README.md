# data/

This directory holds **derived data only**. No large FITS dataset is committed to the repository, and the application
never ships with data it did not download or generate.

```
data/
├── real/        verified real SPHEREx cutouts retrieved from IRSA + manifest.json
├── demo/        synthetic DEMONSTRATION DATA products + manifest.json
├── uploads/     user-supplied local FITS files
├── cache/fits/  HTTP cutout cache keyed by URL hash (safe to delete)
└── README.md
```

## real/ — genuine SPHEREx observations

Produced by:

```bash
python3 scripts/build_demo_bundle.py          # retrieves the verified pair from IRSA
python3 scripts/build_demo_bundle.py --check  # verifies checksums only
```

`real/manifest.json` records, for each cached file: the TAP query executed, the exact download URL, file size,
SHA-256, retrieval timestamp, and the **measured** WCS-WAVE central wavelength and bandwidth, pixel scale, cutout centre,
severe-flag fraction and extension list. If the archive no longer returns the verified observation, the script fails with
an explicit message instead of substituting a different one.

Current contents (~7.7 MB, two files):

| role | observation | date (UTC) | file | λ (µm) | bandwidth (µm) |
|---|---|---|---|---|---|
| epoch 1 | `2026W20_2A_0464_1` | 2026-05-17T10:53:45 | `495d8688953ba51f80530bdb.fits` | 1.42020 | 0.03422 |
| epoch 2 | `2026W26_1A_0178_1` | 2026-06-22T20:21:33 | `e171474e2ab87bfe8c522a3e.fits` | 1.42058 | 0.03423 |

Field: M101, RA 210.80227°, Dec +54.34895°, SPHEREx-D2, 59 × 59 px cutouts at 6.14″/px.

To use a different region, query the archive from **Sky Explorer**; downloads land in `data/cache/fits/` automatically.

## demo/ — labelled synthetic products

Produced by:

```bash
python3 scripts/generate_demo_data.py
```

Two 96 × 96 MEF files with `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, `EPSF` and a TAN WCS. Their headers contain
`DATAMODE = 'DEMONSTRATION DATA - SYNTHETIC'`, and `demo/manifest.json` states the truth model:

- a source that shifts 2.5 px between epochs,
- a source that brightens by 60 %,
- a source that disappears,
- a 3-pixel cosmic-ray-like streak present only in epoch 2,
- a NaN/NONFUNC 2 × 2 patch present in both epochs.

These files exist so the pipeline, interface and tests can be exercised without internet access. They are **not**
SPHEREx observations and must never be presented as such. The UI shows a red `DEMONSTRATION DATA — SYNTHETIC` banner
whenever they are loaded.

## uploads/ — local FITS files

Files uploaded through the interface are stored here with a random prefix and a sanitised name. Their provenance cannot
be verified, so they are labelled `LOCAL FITS UPLOAD · UNVERIFIED PROVENANCE`, and they can only be compared with other
local uploads.

## cache/fits/

Raw IRSA cutout responses, keyed by the SHA-256 of the request URL, with a JSON sidecar containing the URL, checksum,
byte count, content type and retrieval time. Deleting this directory is safe; files are re-downloaded on demand.

## What is deliberately not stored

- Full-frame SPHEREx products (~65–72 MB each). Only cutouts are requested.
- Any data from the archive that failed validation.
- Any result that cannot be traced to two input files with recorded checksums.

## Obtaining SPHEREx data yourself

1. IRSA SPHEREx data access: <https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/>
2. IRSA cutout tutorial (TAP discovery + cutout URLs): <https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/>
3. Browse the archive tree: <https://irsa.ipac.caltech.edu/ibe/data/spherex/>
4. Query interactively with TAP: <https://irsa.ipac.caltech.edu/TAP>

Please cite NASA SPHEREx and IRSA in any downstream use.
