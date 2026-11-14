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
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import config  # noqa: E402
from modules import data_loader, epsf_processor, irsa_client, multi_epoch, pipeline, wcs_alignment  # noqa: E402

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

            print("5b. Fitting one candidate with each real QR3 ePSF model ...")
            if not all("EPSF" in {str(name).upper() for name in product.extension_names} for product in products[:2]):
                failures.append("the verified QR3 reference pair did not both contain the documented EPSF extension")
            elif not outcome.get("candidates"):
                failures.append("the real two-epoch run produced no candidate position for an ePSF fit")
            else:
                epsf_candidate = next(
                    (
                        candidate for candidate in outcome["candidates"]
                        if candidate.quality.get("matched_epoch1") and candidate.quality.get("matched_epoch2")
                    ),
                    outcome["candidates"][0],
                )
                epsf_results = []
                from astropy.coordinates import SkyCoord

                for epoch_index, product in enumerate(products[:2], start=1):
                    epoch = epsf_candidate.metadata.get(f"epoch{epoch_index}", {})
                    try:
                        sky = SkyCoord(float(epoch["ra_deg"]), float(epoch["dec_deg"]), unit="deg", frame="icrs")
                        x_native, y_native = product.wcs.world_to_pixel(sky)
                        fit = epsf_processor.fit_source(product, float(x_native), float(y_native))
                    except Exception as error:
                        fit = {"status": "unavailable", "detail": str(error)}
                    epsf_results.append(fit)
                    if fit.get("status") != "ok":
                        failures.append(f"real epoch {epoch_index} ePSF fit failed: {fit.get('detail', fit.get('status'))}")
                    else:
                        print(
                            f"   epoch {epoch_index} {epoch.get('id')} | candidate {epsf_candidate.candidate_id} | "
                            f"fit {fit['flux_fit']:.6g} {fit['flux_unit']} | "
                            f"reduced chi-square {fit['reduced_chi_square']:.3f} | "
                            f"ePSF zone {fit['epsf_zone_id']} | checksum {fit['product_sha256']}"
                        )
                if len(epsf_results) == 2 and all(row.get("status") == "ok" for row in epsf_results):
                    print("   ePSF fit is supplemental photometry; no classification or sensitivity claim is derived from it")

    print("6. Loading independent dates for a real three-epoch analysis ...")
    by_date: dict[str, dict] = {}
    for observation in sorted(observations, key=lambda row: (row.get("mjd") is None, row.get("mjd") or 0.0)):
        date_utc = str(observation.get("date_utc") or "")
        if not date_utc or not observation.get("cutout_url"):
            continue
        by_date.setdefault(date_utc[:10], observation)
    independent_records = list(by_date.values())[:10]
    loaded_records: list[dict] = []
    loaded_products: dict[str, object] = {}
    retrieval_issues: list[str] = []
    selected_triple: tuple[list[dict], dict[str, object]] | None = None
    multi_cfg = dict(config.PROCESSING_DEFAULTS)
    multi_cfg.update(multi_epoch.DEFAULT_MULTI_EPOCH_CONFIG)
    for record in independent_records:
        try:
            product = data_loader.load_observation_product(record)
        except Exception as error:
            retrieval_issues.append(f"{record.get('obs_id')}: {error}")
            continue
        loaded_records.append(record)
        loaded_products[str(record["id"])] = product
        if len(loaded_records) < 3:
            continue
        matrix = multi_epoch.compatibility_matrix(loaded_records, loaded_products, multi_cfg)
        by_pair = {
            frozenset((str(row["observation_id_1"]), str(row["observation_id_2"]))): row
            for row in matrix
        }
        for triple in combinations(loaded_records, 3):
            ordered = sorted(triple, key=lambda row: (row.get("mjd") is None, row.get("mjd") or 0.0))
            times = [float(row["mjd"]) for row in ordered if row.get("mjd") is not None]
            if len(times) != 3 or min(times[1] - times[0], times[2] - times[1]) < 1.0:
                continue
            pair_checks = [
                by_pair.get(frozenset((str(left["id"]), str(right["id"]))))
                for left, right in combinations(ordered, 2)
            ]
            # The acceptance example is deliberately stricter than a positional
            # association: all three epochs must share a detector, compatible
            # WCS-WAVE sampling, temporal separation, and spatial coverage.
            if all(row and row["photometric_compatible"] for row in pair_checks):
                selected_triple = (
                    ordered,
                    {str(row["id"]): loaded_products[str(row["id"])] for row in ordered},
                )
                break
        if selected_triple:
            break

    if selected_triple is None:
        failures.append(
            "no three independent dates passed every pairwise detector, WCS-WAVE, temporal, and spatial compatibility check "
            f"after loading {len(loaded_records)} product(s)"
        )
    else:
        triple_records, triple_products = selected_triple
        for record in triple_records:
            product = triple_products[str(record["id"])]
            print(f"   {record['date_utc']} | {record.get('obs_id')} | {record.get('band')} | sha256 {product.sha256}")
        analysis = multi_epoch.analyse_multi_epoch(triple_records, triple_products, multi_cfg)
        summary = analysis.get("summary", {})
        if analysis.get("status") != "ok" or summary.get("observations") != 3:
            failures.append(f"three-epoch analysis did not complete: {analysis.get('error', analysis.get('status'))}")
        elif summary.get("detections", 0) == 0:
            failures.append("three-epoch analysis completed but detected no sources in the selected real products")
        else:
            print(
                f"   run {analysis['run_id']} | compatibility pairs {len(analysis['compatibility_matrix'])} | "
                f"detections {summary['detections']} | tracks {summary['tracks']} | "
                f"potential moving-source leads {summary['potential_moving_source_tracks']}"
            )
            print("   motion remains an investigative lead; the script makes no object or discovery claim")
            print("   all selected epoch pairs passed detector, WCS-WAVE, temporal, and spatial compatibility")

    if retrieval_issues:
        print(f"   {len(retrieval_issues)} additional archive product(s) could not be loaded during epoch selection")
        for issue in retrieval_issues:
            print(f"   - {issue}")

    print("\nRESULT:", "PASS" if not failures else "FAIL")
    for failure in failures:
        print("  -", failure)
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
