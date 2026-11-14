"""Live integration check against the real NASA/IPAC IRSA services.

This script is deliberately separate from the pytest suite so that the unit and
integration tests remain runnable offline. It performs the full real workflow:

  1. verify both IRSA endpoints respond,
  2. run the documented TAP discovery query for the verified M101 field,
  3. download two real cutouts (the verified wavelength-matched pair when present),
  4. validate FITS/WCS/WCS-WAVE and run the difference + detection pipeline,
  5. print a summary and exit non-zero on any failure.

Usage:
    python3 scripts/live_archive_check.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config  # noqa: E402
from modules import data_loader, irsa_client, pipeline, wcs_alignment  # noqa: E402

VERIFIED_IDS = ("2026W20_2A_0464_1", "2026W26_1A_0178_1")


def main() -> int:
    failures: list[str] = []

    print("1. Checking IRSA endpoints …")
    status = irsa_client.check_archive_status()
    print(f"   status={status['status']} ({status['detail']}, {status['latency_seconds']} s)")
    for service, detail in status["services"].items():
        print(f"   - {service}: {detail['status']}")
        if detail["status"] != "ok":
            failures.append(f"{service} endpoint failed: {detail.get('error')}")

    print("2. Running the documented TAP discovery query …")
    query = dict(config.DEFAULT_QUERY)
    result = irsa_client.query_archive(query)
    observations = result["observations"]
    print(f"   {len(observations)} real rows · SIA2 coverage {result['coverage_status']['status']}")
    if not observations:
        failures.append("TAP returned no observations for the verified field")

    print("3. Selecting the verified wavelength-matched pair …")
    pair = [obs for obs in observations if obs.get("obs_id") in VERIFIED_IDS]
    if len(pair) < 2:
        print("   verified pair not present in this response; using the first two same-band rows instead")
        pair = observations[:2]
    if len(pair) < 2:
        failures.append("fewer than two observations available")
        print("FAILED:", *failures, sep="\n  - ")
        return 1
    for obs in pair:
        print(f"   {obs['obs_id']} · {obs['date_utc']} · {obs['band']} · release {obs['release']}")

    print("4. Retrieving and validating cutouts …")
    products = []
    for obs in pair:
        try:
            product = data_loader.load_observation_product(obs)
        except Exception as error:
            failures.append(f"cutout retrieval failed for {obs['obs_id']}: {error}")
            continue
        products.append(product)
        wavelength = f"{product.center_wavelength_um:.5f} µm" if product.center_wavelength_um else "unavailable"
        print(
            f"   {product.shape[0]}×{product.shape[1]} px · {product.bunit} · λ={wavelength} · "
            f"{product.file_size_bytes:,} bytes · sha256 {product.sha256[:12]}…"
        )
        if any(check["status"] == "fail" for check in product.validation):
            failures.append(f"FITS validation failed for {obs['obs_id']}")
    if len(products) < 2:
        print("FAILED:", *failures, sep="\n  - ")
        return 1

    print("5. Checking compatibility and running the pipeline …")
    records = pair[:2]
    check = wcs_alignment.check_observation_compatibility(
        records[0], records[1], products[0], products[1], dict(config.PROCESSING_DEFAULTS)
    )
    if not check["compatible"]:
        print("   INCOMPATIBLE:", " | ".join(check["reasons"]))
        failures.append("verified pair is not compatible under the current configuration")
    else:
        print(
            f"   Δλ={check['wavelength_delta_um']:.5f} µm (allowed ≤ {check['wavelength_tolerance_um']:.5f}) · "
            f"Δt={check['time_separation_days']:.2f} d"
        )
        outcome = pipeline.run_pipeline(records[0], records[1], pipeline.processing_config(), "2MASS PSC")
        if outcome["status"] != "ok":
            failures.append(f"pipeline failed: {outcome.get('error')}")
        else:
            print(f"   run {outcome['run_id']} · overlap {outcome['aligned'].overlap_fraction:.1%} · "
                  f"units {outcome['units']} · noise {outcome['difference'].global_noise:.4g}")
            print(f"   candidates: {outcome['candidate_summary']['total']} {outcome['candidate_summary']['counts']}")
            print("   (screening records only — this application makes no discovery claim)")

    print("\nRESULT:", "PASS" if not failures else "FAIL")
    for failure in failures:
        print("  -", failure)
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
