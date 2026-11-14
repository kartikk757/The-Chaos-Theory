# FRONTEND_BACKEND_MAPPING.md

Audit of the existing THE CHAOS THEORY frontend and the backend capability each component requires.

- **Audit date:** 2026-11-13
- **Frontend framework:** Streamlit 1.67 (multi-page, `app.py` + `pages/01…09`), custom CSS design system in `assets/styles.css`, Plotly for scientific charts, Matplotlib for cutout PNGs.
- **Language:** Python 3.13 (3.11+ compatible).
- **Existing backend:** `modules/` (archive client, FITS processor, alignment, normalization, difference, detection, candidate analysis, catalogue matching, visualization, pipeline, state, provenance, models) — 13 modules, ~4.1k lines.
- **New backend layer:** `backend/` — service façade, multi-epoch moving-source engine, scientific report contract, results store. **No frontend file is redesigned.** Existing pages keep their layout, palette, typography, navigation, terminology and interactions.

> **Rule applied:** the backend was shaped to fit the frontend. Where a frontend element required something the science cannot support, the limitation is stated below and the smallest honest modification is proposed instead of a redesign.

---

## 1. Global shell

| Component | File | Current functionality | Required backend functionality | Input | Output | Loading state | Error state | Data visualisation |
|---|---|---|---|---|---|---|---|---|
| Landing hero + CTAs | `app.py` | Identity, mission framing, 2 CTAs, 3 page links | None (static narrative). Attribution text must match verified archive facts | — | — | n/a | n/a | Decorative star field only |
| Persistent sidebar nav | `components/sidebar.py` + `config.NAV_ITEMS` | 9 page links, reset, instrument state block | `backend.observation_manager.inventory()` for counts | session state | counts | Sidebar renders immediately, counts show `0` | Silent; counts stay 0 | Text only |
| Top status bar | `components/navbar.py` | Data mode, IRSA link, processing state, session id | `backend.frontend_adapter.health()` (archive + store reachability) | — | status dicts | `Not checked` until probed | `degraded`/`unavailable` shown, never hidden | Status dots |
| Data-mode banner | `components/status_panel.py` | Distinguishes LIVE / CACHED_REAL / LOCAL / DEMO / EMPTY | `backend.models.observation.data_mode` preserved end-to-end | mode key | banner HTML | n/a | Unknown mode renders `UNKNOWN DATA MODE` | Colour-coded banner |

## 2. Page 01 — Mission Control

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| 5 metric cards | Real session state or `Awaiting data` | `backend.frontend_adapter.summary()` | session state | region, observation count, band, candidate count, processing state | `Awaiting data` empty state | Empty state, never fake values | Glass metric cards |
| Workflow stage strip | 5 stages, done/not-done | Same summary | — | booleans | — | — | Pipeline chips |
| One-click verified pair / demo load | Loads 2 real cached cutouts or synthetic pair | `backend.frontend_adapter.download_observations()` | record ids | validated products | Spinner | Traceback expander | — |
| Session log table | Event log | `backend.utils.logging_config` + `modules.provenance` | — | rows | — | — | Table |

## 3. Page 02 — Sky Explorer

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| RA/Dec/radius/size/date/band/release/survey form | Parses decimal + sexagesimal, validates ranges | `backend.utils.validation.validate_query`, `backend.data_acquisition.search()` | user values | validated query | — | Inline exception + traceback expander | — |
| **QUERY IRSA ARCHIVE** | Live TAP + SIA2 | `backend.data_acquisition.search()` → TAP join, SIA2 footprints, collection discovery | query | observation records + coverage polygons | Spinner | `ArchiveError` surfaced with HTTP status | Sky footprint chart + observation timeline |
| **CHECK ARCHIVE CONNECTION** | Probes TAP + SIA2 | `backend.data_acquisition.health()` | — | status/latency | — | Failure listed per service | — |
| **LOAD VERIFIED CACHED PAIR** | Reads `data/real/manifest.json`, verifies SHA-256 | `backend.frontend_adapter.load_cached_bundle()` | manifest | 2 validated products | Spinner | Checksum mismatch reported | — |
| **LOAD DEMONSTRATION DATA** | Synthetic pair, labelled | `backend.frontend_adapter.load_demo_bundle()` | manifest | 2 synthetic products | — | Missing manifest error | Red DEMO banner |
| Local FITS upload | 1–2 files → validated records | `backend.frontend_adapter.upload_local_fits()` | bytes | LOCAL records | Spinner | Rejected payload explained | — |
| Coverage + timeline charts | Plotly | `modules.visualization` (unchanged) | polygons, records | figures | — | Empty state | Footprint polygons, epoch timeline |
| Query provenance expander | SQL + transport metadata | `backend.data_acquisition.last_query_diagnostics()` | — | SQL text, JSON | — | "no query yet" | Code block |

## 4. Page 03 — Observation Explorer

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| Band/release/dated/cached filters | Client-side filtering | none (pure presentation) | records | subset | — | — | Table |
| Observation table | Real archive metadata | `backend.models.observation.ObservationRecord.to_contract()` | records | table rows | — | `Awaiting data` empty state | Table |
| **LOAD FITS & CHECK COMPATIBILITY** | Downloads, validates, compatibility gate | `backend.frontend_adapter.download_observations()` + `check_compatibility()` | 2 record ids | products + check dict | Spinner per epoch | Retrieval/validation failure with traceback | Validation + compatibility tables |
| Compatibility table | 13 checks incl. WCS-WAVE Δλ | `backend.image_processing.compatibility_report()` | records, products | checks, reasons | — | Refusal reasons shown | Table + warnings |

## 5. Page 04 — Cosmic Image Lab

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| Colormap/contrast/invert controls + reset | Display only | none | — | — | — | — | Rendered arrays |
| Single-epoch viewer + pixel inspector | IMAGE/ZODI/VARIANCE/FLAGS readout at a pixel | `backend.fits_processor.pixel_report()` | product, x, y | value, sky, flags | — | `WCS unavailable` caption | Heatmap |
| Side-by-side / blink / overlay / difference | Plotly panels, fragment timer | `backend.image_processing.registered_pair()` | 2 products | aligned arrays | Fragment timer | Incompatible pair blocked upstream | Heatmaps, blink |
| **RUN ALIGNMENT · NORMALIZATION · DIFFERENCE · DETECTION** | Full two-epoch pipeline | `backend.frontend_adapter.generate_difference_image()` | records, cfg | difference, SNR, candidates | Spinner | Stage table + error | Difference + SNR maps |
| Downloads (FITS, CSV) + provenance JSON | Exact bytes | `modules.provenance` | products | files | — | Missing file hidden | — |

## 6. Page 05 — Detection Engine

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| 4 summary cards | Counts by class | `backend.scientific_analysis.class_summary()` | candidates | counts | `Awaiting` | — | Metric cards |
| Stage table | 9 stages with status | `backend.services.processing_service` result | — | stage rows | — | Failed stage detail | Table |
| Residual map + candidate overlay | Plotly heatmap + markers | `modules.visualization` | arrays, candidates | figure | — | "no residuals" null result | Heatmap + markers |
| Candidate table + CSV | Selection + export | `modules.pipeline.candidates_to_csv` | candidates | CSV | — | Empty state | Table |
| Provenance download | JSON | `modules.provenance` | run | file | — | — | — |

## 7. Page 06 — Candidate Investigation

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| Candidate selector | Dropdown of records | `backend.scientific_analysis.report(candidate)` | id | report dict | — | Empty state | — |
| Position/brightness/displacement cards | Real measurements | report contract (Phase 9 schema) | report | values | — | `not measurable` | Metric cards |
| Candidate + difference cutout panels | Plotly | `modules.visualization` | arrays | figures | — | — | Heatmap, marker |
| "Why was this candidate flagged?" | Generated from measurements | `report["evidence"]` | report | bullet list | — | — | Expander |
| Catalogue cross-match block | 2MASS/AllWISE | `modules.catalogue_matching` | candidate | match dict | — | Service failure shown | Table |
| CSV export | All candidates | `backend.frontend_adapter.export_candidates_csv()` | candidates | CSV | — | — | — |

## 8. Page 07 — Methodology / Page 09 — About

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| Static pipeline explanation, units table, limitations | Narrative | Must stay truthful to `backend/` behaviour | — | — | — | — | Text + table |
| Version table, documentation links | `modules.provenance.software_versions()` | — | — | rows | — | — | Table |

## 9. Page 08 — Data Explorer

| Component | Current functionality | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|---|
| Raw TAP rows, normalised records | Verified metadata only | `backend.models.observation` | records | JSON/tables | — | Empty state | Tables |
| FITS header tabs, FLAGS definitions, validation | Header-level truth | `backend.fits_processor` | product | header dict, flag table | — | "no FLAGS" note | Tables |
| Provenance + directory inventory | Real files and hashes | `backend.observation_manager.inventory()` | — | rows | — | — | Table |

## 10. New page — 10 · Moving-Source Investigation (the only frontend addition)

Batch 2 of this build needed a home for the multi-epoch capability the competition asks for. The existing 9 pages and their visual language are untouched; the new page reuses the same components (`bootstrap_page`, `metric_cards`, Plotly theme, tables, `page_link` footer).

| Component | Required backend | Input | Output | Loading | Error | Visualisation |
|---|---|---|---|---|---|---|
| Epoch set loader (real cached multi-epoch bundle / 3-epoch synthetic) | `backend.frontend_adapter.load_multiepoch_bundle()` | bundle id | ≥2 validated observations | Spinner | Missing/invalid bundle explained | — |
| Compatibility matrix (astrometric vs photometric) | `backend.image_processing.compatibility_matrix()` | epochs | matrix + verdicts | Spinner | Refused pairs listed | Table |
| **RUN MULTI-EPOCH INVESTIGATION** | `backend.frontend_adapter.detect_moving_candidates()` | epochs, config | tracks, reports, stage log | Spinner | Stage table | — |
| Track table | `backend.scientific_analysis.track_table()` | tracks | rows (class, motion, SNR, epochs) | — | Null result message | Table |
| Track detail | `backend.scientific_analysis.report()` | track id | full contract report | — | `Insufficient evidence` label | Metric cards + evidence list |
| Track overlay on reference epoch | `backend.visualization_helpers` (thin) | ref image, positions | figure | — | — | Heatmap + per-epoch markers + arrow |
| Motion diagram (Δξ/Δη vs time) | same | track | scatter + fitted line | — | — | Plotly scatter |
| Save report / load saved runs | `backend.results_store` (local JSON always, Supabase optional) | report | run id | Spinner | "store unavailable" banner | Table |

## 11. Frontend elements that use placeholder or simulated functionality

| Element | Status before this batch | Resolution |
|---|---|---|
| Demo pair (`data/demo/*.fits`) | Synthetic, clearly labelled | Unchanged; still labelled DEMONSTRATION DATA |
| Any candidate without a track | Single pair only | New multi-epoch engine supersedes it for motion; the 2-epoch view remains valid for variability |
| Motion claims with 2 epochs | Displacement only, tolerance-gated | Now requires ≥3 epochs for a *track*; 2-epoch motions get the `2-EPOCH-ONLY` quality flag |
| "Confidence" | Not present (deliberately) | `ranking_score` added as a **documented heuristic**, explicitly not a probability |

## 12. Frontend→backend call map (stable service functions)

| Service function | Frontend caller | Backing implementation |
|---|---|---|
| `search_sky_region(ra, dec, radius, …)` | `pages/02` | `backend.data_acquisition.search()` |
| `download_observations(ids)` | `pages/01`, `pages/03` | `backend.data_acquisition.fetch()` + `backend.fits_processor.load()` |
| `get_observation_metadata(id)` | `pages/03`, `pages/08` | `backend.models.observation` |
| `process_observations(ids)` | `pages/04` | `backend.services.processing_service.run_pair()` |
| `generate_difference_image(a, b)` | `pages/04`, `pages/05` | `backend.difference_engine.` + `modules.difference_engine` |
| `detect_moving_candidates(ids)` | `pages/10` | `backend.candidate_detection.investigate()` |
| `get_candidate_details(candidate_id)` | `pages/06`, `pages/10` | `backend.scientific_analysis.report()` |
| `get_pipeline_status(job_id)` | `pages/04`, `pages/10` | `backend.frontend_adapter.get_pipeline_status()` |
| `health()` / `summary()` | shell, `pages/01` | `backend.frontend_adapter` |
| `save_report()` / `list_saved_reports()` | `pages/10` | `backend.results_store` |

## 13. Limitations where the frontend asks for more than the science can give

1. **"Is this Planet X?"** — cannot be answered from SPHEREx photometry alone. The UI says so, labels results as investigative leads, and never uses a discovery label. The smallest honest modification (already applied) is the explicit disclaimer block on the new page.
2. **Continuous tracks** — SPHEREx revisits a field at intervals of days to months, not continuously. Angular motion is therefore fitted across sparse epochs and quoted with a uncertainty and a χ² consistency value; it is not a continuous trajectory.
3. **Sub-pixel astrometry** — peak positions are quantised at ~6.1″/px. Motion below ~1 px is reported as consistent with zero and flagged `ASTROMETRY-LIMITED`.
4. **Photometric comparison across LVF wavelengths** — refused (exact Δλ gate). The new page distinguishes *astrometric* compatibility (motion) from *photometric* compatibility (brightness) and states which one each result relies on.
5. **Synchronised pan/zoom** between viewer panes — not implemented; the blink comparator is provided (documented limitation, unchanged from batch 1).
