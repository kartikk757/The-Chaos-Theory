# Frontend and backend map

Audit date: 2026-11-14
Frontend: Streamlit multi-page application in `app.py`, `pages/`, and `components/`.
Scientific backend: Python modules in `modules/`. There is no separate backend service layer.

This map describes current executable code. Optional services are user-triggered and can be unavailable without disabling the core IRSA and two-observation workflow.

## Page map

| Page | User-facing capability | Implementing code and data contract |
|---|---|---|
| 01 Mission Control | Session state, data mode, workflow entry points | `components/sidebar.py`, `components/navbar.py`, `modules/state.py` |
| 02 Sky Explorer | IRSA TAP/SIA2 sky search, observation footprints, optional SPHEREx HiPS context and FITS quick-look | `modules/irsa_client.py`, `modules/hips_viewer.py`, `modules/mosaic_cube.py` |
| 03 Observation Explorer | Observation list, two-record comparison selection, separate 2–20 epoch set, FITS metadata, coordinate provenance, pair compatibility | `modules/data_loader.py`, `modules/fits_processor.py`, `modules/wcs_alignment.py`, `modules/coordinates.py` |
| 04 Cosmic Image Lab | Single-image inspection, time-ordered epoch sequence, blink/opacity view, explicit pair registration and difference, uploaded mosaic/HiPS cube spectrum | `modules/visualization.py`, `modules/pipeline.py`, `modules/mosaic_cube.py` |
| 05 Detection Engine | Two-epoch residual pipeline, multi-epoch compatibility matrix and association, optional Source Discovery Tool/SExtractor comparison, controlled injection/recovery report | `modules/pipeline.py`, `modules/multi_epoch.py`, `modules/spherex_sdt_adapter.py`, `modules/validation.py` |
| 06 Candidate Investigation | Candidate measurements, IRSA catalogues, Gaia DR3, user-selected JPL known-object check, ePSF zone diagnostics, candidate-centered cube spectrum | `modules/catalogue_matching.py`, `modules/gaia_client.py`, `modules/jpl_crossmatch.py`, `modules/epsf_processor.py`, `modules/mosaic_cube.py` |
| 07 Methodology | Methods, assumptions, interpretation boundaries | Static explanation; update it whenever the executable analysis changes |
| 08 Data Explorer | Observation/product headers, quality flags, validation, run/external-service provenance, archive query history and resource links | `modules/fits_processor.py`, `modules/provenance.py`, `modules/state.py` |
| 09 About | Application identity and attribution | `config.py`, `modules/provenance.py` |

## Scientific data flow

```text
Sky Explorer
  -> IRSA TAP/SIA2 records
  -> Observation Explorer selects pair and epoch set
  -> data_loader + fits_processor validate products and preserve hashes/WCS
  -> pair: wcs_alignment -> normalization -> difference -> source detection -> catalogue matches
  -> epochs: pairwise compatibility -> peak detections -> sky-coordinate associations -> linear-motion screening
  -> Candidate Investigation shows measured evidence and optional catalogue hypotheses
  -> provenance JSON and candidate CSV remain downloadable
```

The pair pipeline and multi-epoch screening are separate analyses. Spectral quick-look measurements from a mosaic/HiPS cube are also separate from time-domain evidence.

## Service behavior and output boundaries

| Integration | Current behavior | Failure or scientific boundary |
|---|---|---|
| IRSA archive | TAP/SIA2 search, cutout retrieval, FITS checksums and local cache | HTTP failures appear in the UI; there is no fabricated fallback row |
| Gaia DR3 | TAP cone query; Gaia reference-epoch positions propagated to measured candidate epochs; match radius combines catalogue and candidate uncertainty | Linear proper motion only; parallax and radial-velocity perspective effects are not applied; no match is inconclusive |
| ePSF | Selects and plots the documented EPSF spatial zone using WCS and parent-detector metadata | Inspection diagnostic only; not PSF photometry, matched filtering, or PSF-matched differencing |
| IRSA Source Discovery Tool | Optional adapter can invoke the tutorial package plus external SExtractor and compare extracted positions to Photutils detections | Optional executable/configuration required; not an alternative primary pipeline; its spectral-extraction workflow is not implemented here |
| 2MASS PSC / AllWISE | IRSA TAP cross-match with uncertainty-aware separations, multiple counterparts, quality fields and AllWISE motion diagnostics | Catalogue positions are not propagated to each SPHEREx epoch; AllWISE motion fields are diagnostic only |
| JPL SBDB / Horizons | User-supplied known object is resolved and checked at candidate observation epochs | Not a blind small-body search; positional agreement tests a hypothesis and is not identification or discovery |
| Mosaic / HiPS cube | Validated uploaded or fetched FITS cube, WCS/table wavelength support, aperture-minus-annulus quick-look, flux/error diagnostics and SHA-256 | Manual mosaic upload; no automated mosaic job submission; HiPS spectrum is visualization context, not a time-separated observation |
| Multi-epoch association | Sorts by UTC/MJD, reports pairwise spectral/spatial/data-mode checks, associates sky positions, fits local tangent-plane linear motion for supported tracks | A motion label requires at least three detections; two epochs always carry `2-EPOCH-ONLY`; uncertainty floor is conservative, not calibrated |
| Injection/recovery | Reproducible controlled image/kernel experiment with seed, checksums, recovery fractions and null-control peaks | Does not estimate application-wide completeness, calibrated false-positive probabilities, or a sensitivity improvement |

## Interpretation rules

- Residuals and tracks are observational screening records, not object identifications or discoveries.
- Two positions define a displacement only; no orbit or fitted track is inferred from two epochs.
- Non-detections are not treated as transients because coverage, wavelength response, depth and crowding vary.
- Heuristic categories and match decisions are not calibrated posterior probabilities.
- Optional service failure remains visible and leaves the core pair pipeline available.

See [RESOURCE_INTEGRATION_AUDIT.md](../RESOURCE_INTEGRATION_AUDIT.md) for the implementation audit and [README.md](README.md) for setup, current validation status and scientific limitations.
