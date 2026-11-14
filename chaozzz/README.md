# THE CHAOS THEORY

**Explore the Universe. Detect the Unexpected.**

A NASA SPHEREx-powered multi-epoch astronomical exploration platform: select a sky region, retrieve genuine SPHEREx
spectral-image cutouts from the NASA/IPAC Infrared Science Archive (IRSA), align and compare them across epochs, generate
difference images, detect statistically significant residuals, and investigate the resulting **candidate records** with
full provenance and explicit limitations.

> This is an independent hackathon prototype. It is **not** an official NASA or Caltech/IPAC product, it makes **no
> discovery claim of any kind**, and it has no Planet X claim. SPHEREx is not a dedicated planet or asteroid tracker.
> A candidate detection is not a discovery.

---

## 1. What actually works

| Capability | Status | Where |
|---|---|---|
| Cinematic landing page, custom observatory design system | Implemented | `app.py`, `assets/styles.css` |
| Persistent navigation, top status bar, data-mode banners | Implemented | `components/` |
| Mission Control dashboard driven by real session state | Implemented | `pages/01_Mission_Control.py` |
| Coordinate validation (decimal + sexagesimal), query validation | Implemented | `modules/irsa_client.py`, `pages/02_Sky_Explorer.py` |
| Live IRSA TAP discovery (`spherex.plane` ⋈ `spherex.artifact`) | Implemented, live-tested | `modules/irsa_client.py` |
| Live SIA2 footprint coverage (`spherex_dr1`, `spherex_qr3`, deep collections) | Implemented, live-tested | `modules/irsa_client.py` |
| Documented cutout download with retry, size limits and SHA-256 | Implemented, live-tested | `modules/data_loader.py` |
| FITS inspection: IMAGE, FLAGS, VARIANCE, ZODI, EPSF, WCS-WAVE | Implemented, live-tested | `modules/fits_processor.py` |
| Dynamic FLAGS-bit parsing (never masks `SOURCE`/`FULLSAMPLE`) | Implemented | `modules/fits_processor.py` |
| Observation compatibility gate (band, mode, units, WCS, WCS-WAVE Δλ, overlap) | Implemented | `modules/wcs_alignment.py` |
| Reprojection onto a common WCS with mask-aware resampling | Implemented, live-tested | `modules/wcs_alignment.py` (reproject) |
| ZODI subtraction, robust tilted background, robust linear scaling | Implemented | `modules/image_normalization.py` |
| Difference image + VARIANCE-propagated SNR | Implemented, live-tested | `modules/difference_engine.py` |
| Residual detection with peak finding and ± dipole merging | Implemented, live-tested | `modules/source_detection.py` |
| Candidate measurement, association, conservative classification | Implemented, live-tested | `modules/candidate_analysis.py` |
| Catalogue cross-match against IRSA TAP tables (2MASS PSC, AllWISE) | Implemented, live-tested | `modules/catalogue_matching.py` |
| Image viewer: colormaps, contrast, pixel inspector, blink, overlay, full-screen | Implemented | `components/image_viewer.py` |
| Candidate table, detail panel, CSV export, provenance JSON | Implemented | `components/candidate_table.py`, `modules/pipeline.py` |
| Local FITS upload path with unverified-provenance labelling | Implemented | `modules/data_loader.py` |
| Labelled synthetic fallback (DEMONSTRATION DATA) | Implemented | `scripts/generate_demo_data.py` |
| Unit + integration test suite (57 tests) | Passing | `tests/` |
| Multi-epoch tracking beyond two epochs | **Not implemented** (documented limitation) | — |
| PSF-matched differencing | **Not implemented** (documented limitation) | — |

---

## 2. Quick start

```bash
cd the-chaos-theory
python3 -m venv .venv && source .venv/bin/activate      # optional but recommended
pip install -r requirements.txt

# optional: pre-cache the verified real SPHEREx cutout pair (needs internet, ~8 MB)
python3 scripts/build_demo_bundle.py

# run the application
streamlit run app.py
# → http://localhost:8501
```

The application was launched and tested with:

```bash
streamlit run app.py --server.address=0.0.0.0 --server.port=8501 --server.headless=true \
  --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false
```

Health check used during development: `curl -s -o /dev/null -w '%{http_code}' http://127.0.0.1:8501/_stcore/health`
(returns `200`), and every page route returns `200`.

Run the tests:

```bash
python3 -m pytest tests -q                      # 57 offline tests (no network)
python3 scripts/build_demo_bundle.py --check    # verify cached real cutout checksums
python3 scripts/live_archive_check.py           # live end-to-end IRSA check (needs internet)
```

`scripts/live_archive_check.py` is the reproducible live integration test: it verifies both IRSA endpoints, runs the
documented TAP discovery query, downloads the verified pair, validates FITS/WCS/WCS-WAVE, runs the pipeline and prints a
PASS/FAIL summary. Last recorded run: **PASS** (2/2 endpoints, 32 rows, 59×59 px cutouts, Δλ = 0.00038 µm, Δt = 36.39 d,
run `CT-C3D5E7AED29D`, 11 candidate records).

---

## 3. Reproducible judge demonstration

### Path A — verified real SPHEREx observations (primary)

1. Open the app → landing page → **INITIATE EXPLORATION**.
2. **Sky Explorer** → either press **QUERY IRSA ARCHIVE** (live archive, internet required) with the pre-filled M101 field,
   or press **LOAD VERIFIED CACHED PAIR** to use the two real cutouts already downloaded from IRSA.
3. **Observation Explorer** → select the two epochs, press **LOAD FITS & CHECK COMPATIBILITY**. The page prints the measured
   WCS-WAVE central wavelengths, Δλ versus the allowed tolerance, FITS validation checks and overlap status.
4. **Cosmic Image Lab** → inspect each epoch, try side-by-side, blink and opacity overlay, then press
   **RUN ALIGNMENT · NORMALIZATION · DIFFERENCE · DETECTION**.
5. **Detection Engine** → inspect the stage table, diagnostics, residual map with candidate positions, and the candidate table.
6. **Candidate Investigation** → select a candidate; see its position, flux measurements, displacement, quality warnings,
   catalogue cross-match result and the expandable *Why was this candidate flagged?* panel. Export the CSV.
7. **Methodology** and **About** → the scientific reasoning, units, limitations and attributions.

The verified pair shipped in `data/real/manifest.json` is:

| | Epoch 1 | Epoch 2 |
|---|---|---|
| Observation | `2026W20_2A_0464_1` | `2026W26_1A_0178_1` |
| Date (UTC) | 2026-05-17T10:53:45 | 2026-06-22T20:21:33 |
| Band | SPHEREx-D2 | SPHEREx-D2 |
| Measured WCS-WAVE λ | 1.42020 µm | 1.42058 µm |
| Measured bandwidth | 0.03422 µm | 0.03423 µm |
| Δλ | 0.000381 µm (≈1.1 % of the bandwidth) | |
| Δt | 36.39 days | |
| Cutout | 59 × 59 px, 6.14″/px, 4,000,320 bytes each | |

With the default configuration this pair yields a **small set of residual candidate records**. Reproducible run
`CT-C3D5E7AED29D` (deterministic over the two file checksums plus the default configuration) produced 11 records:
2 possible variable-source candidates, 6 unmatched residuals, 3 unclassified detections, 0 moving-source candidates.
The pipeline also found 51 2MASS point sources in the 0.08° cone and matched one candidate within 6″; the rest have no
catalogue counterpart at that radius, which is reported as missing information rather than as a new object. These are
screening records produced from two real observations — **not** discoveries, and the interface says so on every page.

### Path B — cached real cutouts, no internet

Press **LOAD VERIFIED CACHED PAIR** on Mission Control or Sky Explorer. The banner switches to
`REAL SPHEREx · CACHED IRSA CUTOUTS` and the checksums are verified against the manifest.

### Path C — labelled synthetic fallback, no internet, no cached files

Press **LOAD DEMONSTRATION DATA**. Every screen then carries the
`DEMONSTRATION DATA — SYNTHETIC` banner, and the FITS files themselves contain
`DATAMODE = 'DEMONSTRATION DATA - SYNTHETIC'`. The synthetic truth model contains a source that shifts 2.5 px, a source
that brightens 60 %, a source that disappears, a single-epoch cosmic-ray-like streak and a NaN/NONFUNC patch. Running the
pipeline on it demonstrates that the classification logic recovers a moving-source candidate, a variable-source candidate
and an unmatched residual — for **software validation only**.

---

## 4. Architecture

```
app.py                      landing page (presentation only)
config.py                   verified endpoints, bands, FLAGS policy, processing defaults, safety limits
pages/01..09                Streamlit pages: orchestration + presentation, no astronomy math
components/                 reusable UI: theme, navbar, sidebar, metric cards, status panel,
                            image viewer, candidate table
modules/
  irsa_client.py            documented TAP / SIA2 access + coordinate and query validation
  data_loader.py            cutout download, cache, uploads, validated FITS product loading
  fits_processor.py         header inspection, FLAGS parsing, WCS + WCS-WAVE evaluation, safety limits
  wcs_alignment.py          compatibility gate + reprojection onto a common grid
  image_normalization.py    ZODI-aware background removal and robust scaling
  difference_engine.py      D = I2 − I1 with propagated uncertainty
  source_detection.py       significance detection and ± dipole merging
  candidate_analysis.py     aperture photometry, epoch association, classification, reasons
  catalogue_matching.py     IRSA TAP catalogue cone search and matching
  visualization.py          Plotly/Matplotlib rendering (no state, no I/O)
  pipeline.py               end-to-end orchestration + provenance + CSV export
  state.py                  session state, event log, data-mode bookkeeping
  provenance.py             SHA-256, deterministic run IDs, JSON/JSONL logging, versions
  models.py                 typed scientific containers (dataclasses)
assets/                     logo.png, background.png (decorative art), styles.css (self-contained)
data/                       real/ cached verified cutouts, demo/ synthetic products, uploads, cache
outputs/                    runs/<run_id>.json, candidates/<run_id>.csv, events.jsonl, application.log
scripts/                    build_demo_bundle.py, generate_demo_data.py
tests/                      pytest suite (unit + integration)
```

**Separation of concerns.** Pages never compute astronomy: they call `modules/*` and render the results. `visualization.py`
contains no I/O. `fits_processor.py` and `wcs_alignment.py` never touch Streamlit state or the network.

**Data flow.** `IRSA TAP/SIA2 → observation records → cutout download (cached) → validated FitsProduct →
compatibility gate → reprojection → normalization → difference → detection → candidate measurement →
catalogue cross-match → provenance record + CSV export`.

---

## 5. Data sources and verified interfaces

Everything below was read from official documentation and confirmed against the live services on 2026-11-13, before
implementation:

| Interface | Endpoint / form | Use |
|---|---|---|
| IRSA TAP (sync) | `https://irsa.ipac.caltech.edu/TAP/sync` with `QUERY`, `FORMAT=csv` | Observation discovery, catalogue cross-match |
| IRSA SIA2 | `https://irsa.ipac.caltech.edu/SIA?COLLECTION=…&POS=circle RA DEC R&RESPONSEFORMAT=CSV` | Footprint coverage (`s_region`, `s_ra`, `s_dec`, `t_min`) |
| IRSA cutout service | `<artifact.uri>?center=RA,DECd&size=DEG` | Small FITS cutouts of large MEF products |
| SIA2 collections | `spherex_dr1`, `spherex_dr1_deep`, `spherex_qr3`, `spherex_qr3_deep` | Searchable releases (QR2 is superseded and not SIA2-searchable) |
| Catalogue tables | `fp_psc` (2MASS PSC), `allwise_p3as_psd` (AllWISE) | Candidate cross-match context |

The exact SQL used by the application is displayed in **Sky Explorer → Query provenance** and stored in
**Data Explorer**. Discovery query shape:

```sql
SELECT TOP n
  'https://irsa.ipac.caltech.edu/' || a.uri || '?center=RA,DECd&size=SIZE' AS uri,
  a.uri AS raw_uri, p.time_bounds_lower, p.time_bounds_upper, p.time_exposure,
  p.energy_bandpassname, p.energy_bounds_lower, p.energy_bounds_upper,
  p.calibrationlevel, p.quality_flag, p.provenance_version, p.productid,
  p.obsid, p.planeid, p.obs_publisher_did, ...
FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid
WHERE 1 = CONTAINS(POINT('ICRS', RA, DEC), p.poly)
  AND p.time_bounds_lower >= :start_mjd AND p.time_bounds_lower < :end_mjd
  AND a.producttype = 'science' AND a.contenttype = 'image/fits'
  [AND p.energy_bandpassname = 'SPHEREx-D2']
ORDER BY p.time_bounds_lower
```

Documentation consulted (also linked in the UI):

- SPHEREx Data Access — <https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/>
- SPHEREx Data Products (MEF extensions, FLAGS bitmap, LVF bands) — <https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/>
- IRSA cutout tutorial — <https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/>
- IRSA SPHEREx introduction (spectral WCS) — <https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/>
- IRSA TAP documentation — <https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html>
- IRSA cutout parameters — <https://irsa.ipac.caltech.edu/ibe/cutouts.html>
- NASA SPHEREx mission — <https://science.nasa.gov/mission/spherex/>

Data volume policy: full-frame products are ~65–72 MB each, so the application always requests **cutouts** (a few MB).
No real or cached SPHEREx FITS data is committed. The two small, clearly labelled synthetic demo FITS files are included
so the demonstration workflow also works in a fresh cloud deployment; see `data/README.md`.

---

## 6. Methodology (implemented behaviour)

1. **Raw SPHEREx data** — Level-2 spectral-image cutouts with `IMAGE` (MJy/sr, no zodiacal subtraction), `FLAGS`,
   `VARIANCE` ((MJy/sr)²), `ZODI`, `EPSF`/`PSF`, `WCS-WAVE`.
2. **Observation validation** — 2-D array, units, extension-shape agreement, celestial WCS construction, pixel scale,
   positive-variance fraction, named FLAGS bits, and WCS-WAVE evaluation **at the geometric centre of the delivered
   cutout** (the LVF central wavelength varies across the focal plane, and the archive record's coordinates can be the
   centre of the parent observation footprint, not of the cutout).
3. **Compatibility gate** — same detector band, same data mode, same BUNIT, celestial WCS present, and
   Δλ ≤ 25 % of the narrower WCS-WAVE bandwidth (configurable). Cutout centre separation must be ≤ 75 % of the smaller
   cutout width. Failure refuses the comparison with a printed reason.
4. **WCS alignment** — both epochs reprojected onto the epoch-1 celestial grid with `reproject.reproject_interp`
   (bilinear for science and variance, nearest-neighbour for masks). Masked/NaN pixels are zero-filled before
   interpolation and re-masked afterwards. Overlap < 55 % or < 10 % usable pixels raises `AlignmentError`.
5. **Background correction** — documented `ZODI` model subtracted per epoch, then a robust tilted plane (or a
   sigma-clipped constant) fitted to unmasked sky pixels. A robust linear sky scaling coefficient is applied only when
   physically plausible (0.2 ≤ gain ≤ 5).
6. **Image comparison** — `D(x, y) = I₂(x, y) − I₁(x, y)` and `SNR = D / sqrt(VAR₁ + VAR₂)`; a global robust difference
   noise replaces the propagated sigma only when VARIANCE planes are missing, and the fallback is reported.
7. **Source detection** — `photutils.find_peaks` on the SNR map and on its inverse. Because SPHEREx pixels (~6.1″)
   undersample the PSF, a box-car peak finder is used instead of a multi-pixel DAO model; opposite-sign neighbours within
   4 px are merged as a dipole so one shifted object is one record.
8. **Candidate measurement** — annulus-subtracted aperture photometry in both epochs (converted to Jy with the pixel
   solid angle), local difference noise, severe-flag fraction in the 3×3 neighbourhood, independent per-epoch source
   detection, spherical displacement and position angle, Δt, alignment tolerance.
9. **Classification** — `Potential moving-source candidate`, `Potential variable-source candidate`, `Unmatched residual`,
   `Likely image artifact (user review required)`, `Unclassified detection`. No planet/asteroid/discovery labels exist in
   the code path; a unit test asserts their absence.
10. **Catalogue cross-match** — IRSA TAP cone search (default 0.08°) and 6″ matching radius; "no match" is reported as
    missing information, not as a new object.
11. **Provenance** — deterministic `CT-XXXXXXXXXXXX` run ID over input SHA-256 checksums + processing configuration +
    application version; JSON record in `outputs/runs/`, CSV in `outputs/candidates/`, JSONL events in
    `outputs/events.jsonl`.

---

## 7. Integrity rules enforced in code

- Missing/invalid data never becomes a result: loaders raise, pages display the error with a traceback expander.
- `FLAGS != 0` is **never** used as a bad-pixel mask. Bit names/masks come from each FITS header; `SOURCE` and
  `FULLSAMPLE` are informational (in a real M101 cutout, `SOURCE` marks ~99 % of pixels).
- Real archive data, cached real data, local uploads and synthetic data carry distinct `data_mode` labels; mixing them is
  a hard compatibility failure.
- Synthetic products carry `DEMONSTRATION DATA - SYNTHETIC` in their FITS headers and are announced by a red banner.
- Every displayed candidate number is traceable to its two input files, their checksums and the configuration used.
- No network failure is hidden, and no observation is silently substituted.
- Units are explicit everywhere: MJy/sr for surface brightness, Jy for fluxes, arcsec for displacement, days for Δt.

---

## 8. Performance and reliability

- Only cutouts are requested (≈4 MB) instead of ~70 MB parents; downloads are cached by URL hash in `data/cache/fits/`.
- Uploads are limited to 200 MB; downloads to 120 MB; local arrays larger than 4 MPixel are centre-cropped to 512×512 and
  the crop is reported (with the central wavelength marked unverified after cropping).
- Archive metadata cached 30 min (`st.cache_data`), SIA2 15 min, FITS products cached by path + mtime + SHA-256,
  processing results cached by serialized inputs + configuration; session state is never cached.
- HTTP retries (3 attempts with backoff) and explicit `ArchiveError` messages, including the HTTP status and response body
  excerpt.
- Rendering uses downsampled percentiles; heavy work (reprojection, photometry, detection) runs inside cached functions
  separate from rendering.

---

## 9. Testing

```
$ python3 -m pytest tests -q
57 passed
```

- `tests/test_data_loader.py` — extension/unit discovery, FLAGS parsing, NaN handling, upload rejection, checksums,
  local-observation labelling.
- `tests/test_alignment.py` — compatibility gate (band, mode, WCS, Δλ tolerance, missing WCS) and reprojection,
  masking and insufficient-overlap refusal.
- `tests/test_difference.py` — background removal, `D = I₂ − I₁` equality, VARIANCE propagation, NaN preservation,
  robust-noise fallback.
- `tests/test_detection.py` — recovery of a known injected residual, noise-only behaviour, threshold monotonicity,
  dipole merging, and the guard that candidate text never contains discovery language.
- `tests/test_pipeline.py` — end-to-end run, deterministic run IDs, provenance file, CSV export, refusal of an
  incompatible pair, honest failure on a missing file.
- `tests/test_irsa_client.py` — coordinate parsing, query validation, record normalisation, SIA2 collection mapping.

Interface-level verification performed during development:

- `streamlit.testing.v1.AppTest` for `app.py` and all nine pages: 0 exceptions in the empty state.
- `AppTest` with a populated session (demo pair *and* the real cached pair, including a completed pipeline run):
  0 exceptions on Mission Control, Cosmic Image Lab, Detection Engine, Candidate Investigation and Data Explorer.
- Live server: root and every page route return HTTP 200; `/_stcore/health` returns 200.

---

## 10. Known limitations and honest gaps

1. **Two epochs only.** Motion is reported as the displacement between matched counterparts in one pair. Track
   consistency across three or more epochs is not implemented.
2. **No PSF matching.** SPHEREx pixels undersample the PSF, so PSF differences leave residuals. Aperture photometry and
   peak detection are used, not optimal PSF photometry.
3. **Displacement quantisation.** Positions come from peak pixels; sub-pixel centroiding and the resulting ~1-pixel
   quantisation mean small motions are not measurable. The alignment tolerance (1.5 px ≈ 9.2″) is stated on every record.
4. **Catalogue matching is shallow.** 2MASS and AllWISE only; no proper-motion, variability or solar-system ephemeris
   matching, so "no match" is genuinely inconclusive.
5. **Zodiacal handling** relies on the product `ZODI` model; no independent background model is fitted.
6. **WCS-only alignment.** No image-based refinement (cross-correlation) is attempted when the WCS is imperfect.
7. **No synchronized pan/zoom** between the two viewer panels; the blink comparator is provided for registration checks.
8. **Single-pair scope.** The app compares two observations at a time; it does not batch-process a field.
9. **Local uploads** are processed but their provenance cannot be verified, and the interface says so.

---

## 11. Attribution and terms

- Mission data: **NASA SPHEREx**. Archive services: **NASA/IPAC Infrared Science Archive (IRSA), Caltech/IPAC**.
- Cutouts are retrieved live or cached verbatim from IRSA; checksums and retrieval timestamps are recorded in
  `data/real/manifest.json`.
- `assets/background.png` is decorative artwork generated for this project, not SPHEREx imagery.
- This repository is an independent hackathon prototype, not an official NASA or Caltech/IPAC product, and publishes no
  astronomical discovery.
