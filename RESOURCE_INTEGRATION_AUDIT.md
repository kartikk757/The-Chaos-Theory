# Resource Integration Audit

Audit date: 2026-11-14  
Scope: the checked-in `chaozzz/` Streamlit application and its deployed source. This document distinguishes implemented code from optional integrations and experimental analyses. A linked service or tutorial is not represented as available unless the app calls it.

## Summary

The application now supports the original IRSA-driven two-observation workflow and a separate multi-epoch screening path. Users can search real SPHEREx archive records, retrieve and validate cutout FITS products, inspect coordinate provenance, compare a compatible pair, or select a 2–20 observation set for pairwise compatibility checks and candidate-track screening. Optional evidence includes Gaia DR3, 2MASS/AllWISE, ePSF zone diagnostics, a SExtractor comparison adapter, JPL known-object checks, mosaic/HiPS cube quick-look spectra, and controlled injection/recovery.

The integrations have deliberately different evidence boundaries. The Gaia and catalogue matches are supporting context; ePSF plots are inspection only; JPL queries require a user-supplied known object; cube spectra are separate from time-domain measurements; and the injection/recovery tool reports results only for the uploaded control image and kernel. The app makes no Planet X, discovery, calibrated probability, or sensitivity-improvement claim.

## Capability status

| Area | Status | Repository evidence and limits |
|---|---|---|
| Existing Streamlit interface and navigation | Implemented | `chaozzz/app.py`, nine pages in `chaozzz/pages/`, and reusable components. No new page or service façade was introduced. |
| IRSA SPHEREx observation discovery | Implemented | `modules/irsa_client.py` performs TAP discovery and SIA2 coverage lookup; failures, TAP diagnostics and a session query history are surfaced. |
| FITS download, validation, checksum and cache | Implemented | `modules/data_loader.py` validates archive URLs, size-limits and hashes downloads, and distinguishes real archive data from local uploads and demo data. |
| Coordinate handling and provenance | Implemented | `modules/coordinates.py`, `data_loader.py`, `fits_processor.py`, and `wcs_alignment.py` distinguish requested field position, archive record position, delivered cutout WCS centre, and wavelength sample. Pair checks use FITS WCS centres. |
| Original two-observation image pipeline | Implemented | `modules/pipeline.py`, `wcs_alignment.py`, `image_normalization.py`, `difference_engine.py`, `source_detection.py`, and `candidate_analysis.py` validate, align, compare, detect and report candidate measurements. |
| Multi-epoch source association | Implemented; experimental | `modules/multi_epoch.py` sorts by UTC/MJD, reports pairwise data-mode/spectral/spatial compatibility, detects peaks, associates positions in a tangent plane, and fits linear motion where supported. Motion classification needs at least three detections. Its positional error floor is not calibrated. |
| Two-epoch interpretation | Explicitly limited | The multi-epoch engine accepts two observations to report pairwise displacement context but attaches `2-EPOCH-ONLY`; it cannot label a fitted motion track or infer an orbit from two positions. |
| Gaia DR3 cross-match | Implemented; optional online service | `modules/gaia_client.py` queries ESA TAP, retrieves astrometry and quality fields, propagates Gaia positions using proper motion to each measured epoch, and uses combined positional uncertainty for matching. Missing proper motion over a longer baseline is reported as unavailable. Parallax and radial-velocity perspective effects are not modelled. |
| ePSF integration | Zone inspection implemented; photometry not implemented | `modules/epsf_processor.py` reads the EPSF extension, maps image coordinates through FITS WCS into the parent detector zone, selects and plots the oversampled model, and reports wavelength/shape diagnostics. Candidate fluxes are not ePSF fits; no PSF-matched differencing or sensitivity gain is claimed. |
| SPHEREx Source Discovery Tool | Optional extraction adapter implemented | `modules/spherex_sdt_adapter.py` checks for the tutorial package, SExtractor and configuration files, runs extraction when available, and compares source positions with Photutils. It is an optional complementary route. The tutorial's separate spectral-extraction workflow is not implemented. |
| 2MASS PSC and AllWISE matching | Implemented; improved | `modules/catalogue_matching.py` uses uncertainty-aware radii, retains multiple counterparts, records catalogue/query and quality information, and labels AllWISE apparent-motion values as diagnostics. Catalogue coordinates are not epoch-propagated in this module. |
| JPL SBDB / Horizons | Implemented for a supplied known-object hypothesis | `modules/jpl_crossmatch.py` resolves a user-entered identifier, retrieves object/orbit details, requests Horizons positions at candidate UTC-derived times and stores angular separations. It does not search blindly for unknown objects. Results and failures are appended to the run provenance JSON. |
| Spectral cube quick-look | Implemented for uploaded and selected HiPS cubes | `modules/mosaic_cube.py` validates 3-D FITS data, spatial/spectral WCS or wavelength table, channel products and units, and calculates aperture-minus-annulus spectra with diagnostics. The UI supports manual mosaic upload and candidate-centered extraction. It does not submit automated mosaic jobs. |
| HiPS sky context | Optional browser and FITS integration | `modules/hips_viewer.py` reads public HiPS properties, calculates NESTED HEALPix indices with optional `hpgeom` or `astropy-healpix`, constructs tile URLs, and bounds tile downloads. Aladin Lite is browser-loaded. HiPS is visualization context and not time-separated evidence. |
| Multi-epoch temporal viewer | Implemented | `pages/04_Cosmic_Image_Lab.py` orders loaded epochs, supports previous/next and timed sequence, opacity and registered difference modes, candidate markers, and time/coordinate labels. Pairwise difference requires the selected pair to pass compatibility. |
| False-positive and injection/recovery evaluation | Controlled experiment implemented; not globally calibrated | `modules/validation.py` evaluates a chosen image and detector-sampled kernel at selected levels, recording seed, hashes, recoveries, null peaks and assumptions. Results are not a global completeness curve or calibrated false-positive probability. |
| Optional dependencies and network handling | Implemented with limits | Core dependencies remain in `requirements.txt`; HEALPix support is optional in `requirements-optional.txt`. HTTP timeouts, bounded retries, caching and process-local rate spacing are applied to IRSA, Gaia, JPL and HiPS calls. SExtractor and its tutorial configuration remain user-installed optional components. |
| Candidate and run provenance | Implemented | Existing pair runs record input IDs, coordinate roles, product hashes, configuration, software versions and external-catalogue diagnostics. Multi-epoch runs receive separate provenance records; a JPL hypothesis check updates the pair-run JSON. |
| Test suite and live archive validation | Validation pending in this environment | New tests cover the integrations and Python compilation succeeds. Pytest cannot import the compiled `astropy.wcs._wcs` extension because Windows Application Control blocks that DLL in this environment. No full suite or live archive pass is claimed. The live check is `chaozzz/scripts/live_archive_check.py`. |

## Required development order

The requested sequence was followed in the implementation:

1. Audit the existing code and document gaps — recorded above.
2. Fix coordinate and provenance consistency — distinct coordinate roles retained and displayed.
3. Build multi-epoch processing — separate epoch-set workflow and explicit `2-EPOCH-ONLY` handling.
4. Integrate Gaia DR3 — epoch-propagated, uncertainty-aware matching.
5. Integrate ePSF processing — selection and diagnostics only; measurement remains evidence-gated.
6. Add the Source Discovery Tool adapter — optional SExtractor path with fallback to Photutils.
7. Improve catalogue matching — uncertainty, multiple matches, provenance and motion diagnostics.
8. Add JPL known-object investigation — SBDB/Horizons at observed candidate epochs.
9. Add mosaic cube analysis — FITS validation and quick-look spectrum extraction.
10. Add HiPS visualization — optional browser viewer and bounded cube retrieval.
11. Improve temporal blinking — epoch sequencing, speed, navigation, overlays and registered comparison.
12. Run complete tests and end-to-end validation — pending until dependencies are available and commands complete.

## Frontend mapping

`chaozzz/FRONTEND_BACKEND_MAPPING.md` now documents the nine actual Streamlit pages and their current module contracts. Previous claims of a separate `backend/` service façade, page 10, Supabase-backed results store, and additional page flows are not treated as implemented functionality.

## Primary scientific references

- [SPHEREx archive and data-access documentation](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/)
- [IRSA SPHEREx cutout tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/)
- [IRSA SPHEREx ePSF tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/)
- [IRSA SPHEREx Source Discovery Tool tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-source-discovery-tool-demo/)
- [IRSA SPHEREx mosaic tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-mosaic/)
- [IRSA SPHEREx HiPS cube tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-hipscube/)
- [Gaia DR3 source data model](https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html)
- [JPL Small-Body Database API](https://ssd-api.jpl.nasa.gov/doc/sbdb.html) and [Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html)
- [Aladin Lite documentation](https://cds-astro.github.io/aladin-lite/index.html)

## Scientific and product constraints

- A candidate or unmatched catalogue result is not an identification or discovery.
- Two detections define a displacement only. At least three detections are required for a motion classification; no orbital trajectory is inferred.
- The multi-epoch astrometric floor is conservative and currently uncalibrated against SPHEREx residuals.
- Gaia proper-motion propagation is linear and does not model parallax or radial-velocity perspective acceleration.
- The ePSF is displayed for diagnostics only. PSF photometry, matched filtering and PSF matching require validated real-data and injection/recovery evidence first.
- JPL Horizons is used only for a user-specified known-object hypothesis.
- A mosaic or HiPS spectrum is a separate quick-look measurement, not time-domain evidence. IRSA's HiPS spectrum remains visualization-only context.
- Injection/recovery outputs describe only the supplied image, kernel and settings; they are not calibrated probabilities or sensitivity claims.
- Optional services must show errors and leave the core pair workflow usable.
