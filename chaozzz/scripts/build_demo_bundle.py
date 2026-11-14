"""Retrieve and cache a verified real SPHEREx cutout pair for offline demonstration.

Run from the project root:

    python3 scripts/build_demo_bundle.py            # retrieve the verified pair
    python3 scripts/build_demo_bundle.py --check    # verify cached files only

Only documented IRSA services are used, and every value written to the manifest
is measured from the retrieved files or taken from the archive response.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config  # noqa: E402
from modules import data_loader, irsa_client, provenance  # noqa: E402

VERIFIED_PAIR = [
    {
        "obs_id": "2026W20_2A_0464_1",
        "role": "epoch 1",
        "expected_date": "2026-05-17",
        "note": "Verified by WCS-WAVE sampling at the M101 cutout centre during development.",
    },
    {
        "obs_id": "2026W26_1A_0178_1",
        "role": "epoch 2",
        "expected_date": "2026-06-22",
        "note": "Wavelength-matched partner of epoch 1 to within ~1% of the WCS-WAVE bandwidth.",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def retrieve() -> dict:
    query = dict(config.DEFAULT_QUERY)
    result = irsa_client.query_archive(query)
    observations = result["observations"]
    if not observations:
        raise SystemExit("IRSA returned no observations for the default query; nothing was cached.")
    selected = []
    for wanted in VERIFIED_PAIR:
        match = next((obs for obs in observations if obs.get("obs_id") == wanted["obs_id"]), None)
        if match is None:
            raise SystemExit(
                f"Verified observation {wanted['obs_id']} is not present in the current archive response. "
                "Not substituting another observation; re-run with a wider date range or investigate the archive."
            )
        selected.append((wanted, match))
    entries = []
    for wanted, record in selected:
        fetched = data_loader.fetch_cutout_path(record["cutout_url"], force=True)
        target = config.REAL_DATA_DIR / Path(fetched["path"]).name
        shutil.copyfile(fetched["path"], target)
        product = data_loader.load_fits_product(target, record, center_ra=record.get("ra_deg"), center_dec=record.get("dec_deg"))
        entries.append({
            **{key: value for key, value in record.items() if key != "metadata"},
            "role": wanted["role"],
            "cached_path": str(target.relative_to(config.ROOT)),
            "derived_sha256": sha256(target),
            "file_size_bytes": target.stat().st_size,
            "download_url": record["cutout_url"],
            "tap_endpoint": config.IRSA_TAP_SYNC,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "measured": {
                "shape": list(product.shape),
                "bunit": product.bunit,
                "center_wavelength_um": product.center_wavelength_um,
                "center_bandwidth_um": product.center_bandwidth_um,
                "pixel_scale_arcsec": product.pixel_scale_arcsec,
                "cutout_center_sky": list(product.center_sky) if product.center_sky else None,
                "wavelength_sample_pixel": product.metadata.get("wavelength_sample_pixel"),
                "severe_flag_fraction": product.metadata.get("flag_diagnostics", {}).get("bad_pixel_fraction"),
                "source_flag_fraction": product.metadata.get("flag_diagnostics", {}).get("source_flag_fraction"),
                "extension_names": product.extension_names,
                "primary_version": product.primary_header.get("VERSION"),
            },
        })
    manifest = {
        "dataset": "THE CHAOS THEORY verified real SPHEREx cutout pair",
        "data_mode": config.DATA_MODES["CACHED_REAL"],
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "query": query,
        "tap_query": result["tap_diagnostics"]["query"],
        "tap_endpoint": config.IRSA_TAP_SYNC,
        "coverage_endpoint": config.IRSA_SIA2,
        "observations": entries,
        "attribution": {
            "mission": "NASA SPHEREx",
            "archive": "NASA/IPAC Infrared Science Archive (IRSA), Caltech/IPAC",
            "note": "Cutouts retrieved from the official IRSA cutout service. This project is not an official NASA product.",
        },
        "integrity_note": "Re-verify checksums with 'python3 scripts/build_demo_bundle.py --check' before a demonstration.",
        "software_versions": provenance.software_versions(),
    }
    provenance.write_json(config.REAL_DATA_DIR / "manifest.json", manifest)
    return manifest


def check() -> int:
    path = config.REAL_DATA_DIR / "manifest.json"
    if not path.exists():
        print("No manifest found. Run 'python3 scripts/build_demo_bundle.py' first.")
        return 1
    manifest = provenance.read_json(path)
    failures = 0
    for entry in manifest["observations"]:
        file_path = ROOT / entry["cached_path"]
        if not file_path.exists():
            print(f"MISSING  {file_path}")
            failures += 1
            continue
        checksum = sha256(file_path)
        status = "OK" if checksum == entry["derived_sha256"] else "CHECKSUM MISMATCH"
        if status != "OK":
            failures += 1
        print(f"{status:16s} {file_path.name}  λ={entry['measured']['center_wavelength_um']:.5f} µm  {file_path.stat().st_size:,} bytes")
    print(f"manifest generated {manifest['generated_utc']} · {len(manifest['observations'])} observations · {failures} problem(s)")
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify cached files instead of retrieving")
    options = parser.parse_args()
    if options.check:
        raise SystemExit(check())
    manifest = retrieve()
    for entry in manifest["observations"]:
        measured = entry["measured"]
        print(
            f"{entry['role']}: {entry['date_utc']} {entry['band']} {entry['obs_id']} "
            f"λ={measured['center_wavelength_um']:.5f} µm Δλ={measured['center_bandwidth_um']:.5f} µm "
            f"→ {entry['cached_path']} ({entry['file_size_bytes']:,} bytes)"
        )
