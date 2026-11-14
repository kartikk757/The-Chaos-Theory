# Judge demonstration script — THE CHAOS THEORY

Target duration: 5–7 minutes. Every step below uses only functionality that exists in this build.

## 0. Before you start (30 s)

```bash
pip install -r requirements.txt
python3 scripts/build_demo_bundle.py --check   # confirms the two real cutouts are intact
streamlit run app.py
```

If `--check` reports missing files, run `python3 scripts/build_demo_bundle.py` with internet access, or use Path C below.

## 1. Landing (30 s)

Point out the hero statement, the SPHEREx/multi-epoch framing, and the statement that candidate detection is not
discovery. Click **INITIATE EXPLORATION**.

## 2. Sky Explorer (60 s)

- The M101 preset is a verified SPHEREx field (RA 210.80227°, Dec +54.34895°).
- Press **QUERY IRSA ARCHIVE** for live metadata, or **LOAD VERIFIED CACHED PAIR** offline.
- Show the requirement/value separation: the requested circle is *not* coverage; cyan polygons are real SIA2 footprints.
- Open *Query provenance* to show the exact SQL and the archive response metadata.

## 3. Observation Explorer (60 s)

- Show the real observation rows with dates, band, release, quality and file availability.
- Select the two epochs and press **LOAD FITS & CHECK COMPATIBILITY**.
- Read the measured Δλ (0.000381 µm) against the allowed tolerance, and the FITS validation table. Emphasise that a pair
  with different wavelengths, bands or data modes is *refused*, not silently compared.

## 4. Cosmic Image Lab (90 s)

- Inspect epoch 1 and epoch 2, then cycle side-by-side → blink → opacity overlay.
- Press **RUN ALIGNMENT · NORMALIZATION · DIFFERENCE · DETECTION**.
- Show the stage table, then the difference image and SNR map. State clearly: a residual is not proof of motion.

## 5. Detection Engine + Candidate Investigation (120 s)

- Stage list, detection diagnostics (peak finding because SPHEREx undersamples the PSF, ± dipoles merged).
- Candidate table → select a candidate.
- Read out: position, aperture fluxes in µJy, displacement in arcsec, Δt, flag fraction, catalogue cross-match status.
- Open **WHY WAS THIS CANDIDATE FLAGGED?** — the text is generated from those measurements.
- Export the CSV; show the run ID and the provenance JSON.

## 6. Methodology and About (45 s)

- The eight-stage diagram, the explanation of why candidates are not discoveries, the units table.
- Attribution, verified documentation links, and the explicit list of known limitations.

## 7. Fallback if the archive is unreachable

- **LOAD VERIFIED CACHED PAIR** — real cutouts, no network, banner reads
  `REAL SPHEREx · CACHED IRSA CUTOUTS`.
- **LOAD DEMONSTRATION DATA** — synthetic; a red `DEMONSTRATION DATA — SYNTHETIC` banner appears on every page and the FITS
  headers say so. Use it to show the classification logic recovering a moving, a variable and an unmatched candidate.
  Say explicitly that these are synthetic test products.

## Questions you may be asked

| Question | Answer |
|---|---|
| Did you discover anything? | No, and the application never claims to. It produces screening records from two real observations. |
| Is this real SPHEREx data? | Yes for Paths A and B — cutouts from IRSA with recorded URLs and SHA-256. Path C is labelled synthetic. |
| Why not use the full 70 MB images? | Cutouts are scientifically sufficient for a small field and fit the hackathon hardware; the app documents the crop. |
| How do you know the images are comparable? | Band, data mode, units, celestial WCS, and a measured WCS-WAVE Δλ tolerance gate; failing pairs are refused. |
| How would you confirm a candidate? | More epochs (track consistency), independent imagery, deeper catalogue/photometric context, and human review. |
