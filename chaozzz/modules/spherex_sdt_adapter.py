"""Optional adapter for the official IRSA SPHEREx Source Discovery Tool.

The core app remains usable without the package or a SExtractor executable.
This adapter uses the IRSA package's packaged SExtractor configuration files
and invokes the installed executable directly, keeping this integration
separate from the Photutils residual detector.
"""
from __future__ import annotations

import importlib.util
import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time

from modules.models import FitsProduct


# IRSA's Source Discovery Tool tutorial documents these DR1 calibration paths
# as of November 2026. The collection is release-specific and must be revisited
# when IRSA publishes a newer solid-angle pixel map collection.
SDT_SAPM_COLLECTION = "cal-sapm-v3-2026-191"
SDT_SAPM_S3 = tuple(
    f"nasa-irsa-spherex/dr1/solid_angle_pixel_map/{SDT_SAPM_COLLECTION}/{detector}/"
    f"solid_angle_pixel_map_D{detector}_spx_{SDT_SAPM_COLLECTION}.fits"
    for detector in range(1, 7)
)
SDT_BAD_FLAGS = (
    "TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "DICHROIC", "MISSING_DATA",
    "HOT", "COLD", "NONLINEAR", "CROSSTALK", "GHOST", "GHOST_FPA", "GHOST_EXT",
    "STREAK", "PERSIST", "BLOOM", "SNOWBALL", "HALO", "SATELLITE_HALO",
)


class SDTUnavailable(RuntimeError):
    """Raised when the optional package or SExtractor binary is unavailable."""


def availability() -> dict[str, Any]:
    package = importlib.util.find_spec("spherex_source_discovery_tool") is not None
    executable = shutil.which("sex") or shutil.which("source-extractor")
    return {
        "status": "available" if package and executable else "unavailable",
        "package_available": package,
        "sextractor_executable": executable,
        "detail": (
            "Official IRSA package and SExtractor executable are available."
            if package and executable
            else "Optional path unavailable: install the IRSA spherex_source_discovery_tool package and SExtractor; the built-in Photutils path remains active."
        ),
        "source": "IRSA SPHEREx Source Discovery Tool tutorial",
    }


def spectrum_availability() -> dict[str, Any]:
    """Report availability of the optional SDT aperture-photometry helpers."""
    try:
        from spherex_source_discovery_tool import aperture_photometry as _aperture_photometry  # noqa: F401
    except Exception as exc:
        return {
            "status": "unavailable",
            "detail": f"Optional IRSA SDT aperture-photometry module is unavailable: {exc}",
            "source": "IRSA SPHEREx Source Discovery Tool tutorial",
        }
    return {
        "status": "available",
        "detail": "IRSA SDT aperture-photometry components are available; calibration maps are checked when extraction runs.",
        "source": "IRSA SPHEREx Source Discovery Tool tutorial",
    }


def _configuration_files() -> dict[str, Path]:
    try:
        from spherex_source_discovery_tool.source_extraction import get_sextractor_file
    except Exception as exc:
        raise SDTUnavailable(f"Could not import the official SPHEREx Source Discovery Tool package: {exc}") from exc
    output = {}
    for filename in ("default_sdt.sex", "default_sdt.param", "default.nnw", "default.conv"):
        try:
            path = Path(get_sextractor_file(filename, package="spherex_source_discovery_tool")).resolve()
        except Exception as exc:
            raise SDTUnavailable(f"Could not locate packaged SExtractor file {filename}: {exc}") from exc
        if not path.is_file():
            raise SDTUnavailable(f"Packaged SExtractor file does not exist: {filename}")
        output[filename] = path
    return output


def _find_column(row: dict, *names: str, default=None):
    lowered = {str(key).lower(): key for key in row}
    for name in names:
        key = lowered.get(name.lower())
        if key is not None:
            return row[key]
    return default


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def extract_sources(product: FitsProduct, timeout_seconds: int = 45) -> dict[str, Any]:
    """Run optional SExtractor on the validated science array and return a table."""
    state = availability()
    if state["status"] != "available":
        return {**state, "status": "unavailable", "sources": [], "count": 0}
    image = np.asarray(product.science, dtype=np.float32)
    if image.ndim != 2 or min(image.shape) < 8 or not np.any(np.isfinite(image)):
        return {"status": "invalid_input", "detail": "A finite two-dimensional validated science image is required.", "sources": [], "count": 0}
    finite = np.isfinite(image)
    fill = float(np.nanmedian(image[finite]))
    image = np.where(finite, image, fill).astype(np.float32)
    try:
        configs = _configuration_files()
        with tempfile.TemporaryDirectory(prefix="chaos_sdt_") as temp_name:
            directory = Path(temp_name)
            image_path = directory / "validated_science.fits"
            catalog_path = directory / "sources.cat"
            # A fresh primary-HDU image prevents writing back to the input file.
            input_header = fits.Header()
            if product.header.get("BUNIT"):
                input_header["BUNIT"] = str(product.header["BUNIT"])
            if product.wcs is not None:
                input_header.extend(product.wcs.to_header(relax=True), update=True)
            fits.PrimaryHDU(data=image, header=input_header).writeto(image_path)
            for name, source in configs.items():
                shutil.copyfile(source, directory / name)
            executable = state["sextractor_executable"]
            command = [
                executable,
                str(image_path),
                "-c", str(directory / "default_sdt.sex"),
                "-CATALOG_NAME", str(catalog_path),
                "-PARAMETERS_NAME", str(directory / "default_sdt.param"),
                "-FILTER_NAME", str(directory / "default.conv"),
                "-STARNNW_NAME", str(directory / "default.nnw"),
                "-CHECKIMAGE_TYPE", "NONE",
            ]
            completed = subprocess.run(
                command,
                cwd=directory,
                capture_output=True,
                text=True,
                timeout=max(1, min(int(timeout_seconds), 120)),
                check=False,
            )
            if completed.returncode != 0 or not catalog_path.is_file():
                return {
                    "status": "failed",
                    "detail": (completed.stderr or completed.stdout or "SExtractor produced no catalogue")[-1200:],
                    "return_code": completed.returncode,
                    "sources": [],
                    "count": 0,
                    "configuration": {key: value.name for key, value in configs.items()},
                }
            table = Table.read(catalog_path, format="ascii.sextractor")
            sources = []
            for row in table:
                item = {str(name): (row[name].item() if hasattr(row[name], "item") else row[name]) for name in table.colnames}
                x = _find_column(item, "X_IMAGE", "XWIN_IMAGE")
                y = _find_column(item, "Y_IMAGE", "YWIN_IMAGE")
                try:
                    x0, y0 = float(x) - 1.0, float(y) - 1.0  # SExtractor pixels are one-based.
                except (TypeError, ValueError):
                    x0 = y0 = None
                if x0 is not None and y0 is not None and product.wcs is not None:
                    try:
                        coordinate = product.wcs.pixel_to_world(x0, y0).icrs
                        item["ra_deg"] = float(coordinate.ra.deg)
                        item["dec_deg"] = float(coordinate.dec.deg)
                    except Exception:
                        item["ra_deg"] = item["dec_deg"] = None
                else:
                    item["ra_deg"] = item["dec_deg"] = None
                item["x_pixel"] = x0
                item["y_pixel"] = y0
                sources.append(item)
            return {
                "status": "ok",
                "sources": sources,
                "count": len(sources),
                "configuration": {
                    "threshold_sigma": None,
                    "configuration_files": {key: str(value) for key, value in configs.items()},
                    "timeout_seconds": int(timeout_seconds),
                    "input": "validated IMAGE minus ZODI array; non-finite pixels replaced by the finite-image median",
                },
                "stdout_tail": completed.stdout[-800:],
                "stderr_tail": completed.stderr[-800:],
                "limitation": "Independent SExtractor extraction; it does not replace the Photutils residual pipeline or constitute validated transient photometry.",
            }
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "detail": f"SExtractor exceeded {timeout_seconds} seconds.", "sources": [], "count": 0}
    except SDTUnavailable as exc:
        return {**availability(), "status": "unavailable", "detail": str(exc), "sources": [], "count": 0}
    except Exception as exc:
        return {"status": "failed", "detail": str(exc), "sources": [], "count": 0}


def compare_with_photutils(sdt_sources: list[dict], photutils_detections: list[dict], match_radius_px: float = 2.0) -> dict[str, Any]:
    """Report one-to-one spatial agreement without ranking either extractor."""
    available = []
    for source in sdt_sources:
        try:
            x, y = float(source["x_pixel"]), float(source["y_pixel"])
            if np.isfinite(x) and np.isfinite(y):
                available.append((x, y))
        except (KeyError, TypeError, ValueError):
            continue
    used: set[int] = set()
    matches = []
    for detection in photutils_detections:
        x = float(detection.get("xcentroid", detection.get("x_peak", np.nan)))
        y = float(detection.get("ycentroid", detection.get("y_peak", np.nan)))
        options = [(float(np.hypot(x - sx, y - sy)), index) for index, (sx, sy) in enumerate(available) if index not in used]
        if not options:
            continue
        separation, index = min(options)
        if separation <= float(match_radius_px):
            used.add(index)
            matches.append({"photutils_index": len(matches), "sdt_index": index, "separation_px": separation})
    return {
        "match_radius_px": float(match_radius_px),
        "photutils_count": len(photutils_detections),
        "sdt_count": len(available),
        "matched_pairs": len(matches),
        "photutils_only_count": max(0, len(photutils_detections) - len(matches)),
        "sdt_only_count": max(0, len(available) - len(matches)),
        "matches": matches,
        "interpretation": "Spatial agreement diagnostic only; these extraction paths differ, and their agreement count is not a calibrated completeness estimate.",
    }


def extract_epoch_spectrum(
    products: list[FitsProduct],
    source_positions: list[tuple[float, float]],
    *,
    source_id: str,
    aperture_radius_px: float = 3.0,
    sapm_s3: tuple[str, ...] | list[str] = SDT_SAPM_S3,
) -> dict[str, Any]:
    """Extract official SDT-style aperture measurements from selected epochs.

    This is intentionally for a small, already validated set of loaded products.
    Each observation uses its own sky coordinate so candidate motion is not
    erased by forcing a single stationary position across epochs. The returned
    wavelength samples remain time-stamped measurements, not a variability-
    corrected or continuous spectrum.
    """
    try:
        aperture_radius_px = float(aperture_radius_px)
    except (TypeError, ValueError, OverflowError):
        return {"status": "invalid_input", "detail": "The aperture radius must be numeric."}
    if not products or len(products) != len(source_positions):
        return {"status": "invalid_input", "detail": "Provide at least one validated product and exactly one ICRS position per product."}
    if not np.isfinite(aperture_radius_px) or aperture_radius_px <= 0 or aperture_radius_px > 5:
        return {"status": "invalid_input", "detail": "The SDT aperture radius must be greater than zero and no larger than 5 pixels."}
    for position in source_positions:
        try:
            ra_deg, dec_deg = float(position[0]), float(position[1])
        except (TypeError, ValueError, IndexError, OverflowError):
            return {"status": "invalid_input", "detail": "Every epoch requires a finite two-value ICRS position."}
        if not np.isfinite(ra_deg) or not np.isfinite(dec_deg) or not 0.0 <= ra_deg < 360.0 or not -90.0 <= dec_deg <= 90.0:
            return {"status": "invalid_input", "detail": "Every epoch requires RA in [0, 360) and Dec in [-90, 90] degrees."}

    try:
        from astropy import units as u
        from astropy.wcs import WCS
        from photutils.aperture import CircularAperture, aperture_photometry
        from spherex_source_discovery_tool.aperture_photometry import (
            extract_dqf,
            extract_photometry,
            fetch_sapm,
            get_AB,
            sa_mjy_scale_map,
        )
    except Exception as exc:
        return {
            "status": "unavailable",
            "detail": f"The optional IRSA SPHEREx Source Discovery Tool photometry components are unavailable: {exc}",
            "measurements": [],
        }

    try:
        # The official helper downloads the six detector maps once; its
        # interpolation routines are then reused for every selected epoch.
        sapm_maps, sapm_interpolators = fetch_sapm(list(sapm_s3))
    except Exception as exc:
        return {
            "status": "unavailable",
            "detail": f"Could not load the official SPHEREx solid-angle pixel maps: {exc}",
            "measurements": [],
            "sapm_collection": SDT_SAPM_COLLECTION if tuple(sapm_s3) == SDT_SAPM_S3 else "caller-supplied",
        }

    output: list[dict[str, Any]] = []
    for product, position in zip(products, source_positions):
        measurement: dict[str, Any] = {
            "observation_id": str(product.product_id),
            "product_sha256": str(product.sha256),
            "source_position_icrs_deg": list(position),
        }
        try:
            if not bool(getattr(product, "is_real_spherex", False)):
                raise ValueError("SDT calibration is enabled only for validated real SPHEREx archive products.")
            ra_deg, dec_deg = float(position[0]), float(position[1])
            if not np.isfinite(ra_deg) or not np.isfinite(dec_deg) or not -90.0 <= dec_deg <= 90.0:
                raise ValueError("A finite ICRS coordinate with declination in [-90, 90] degrees is required.")
            path = Path(product.path)
            if not path.is_file():
                raise ValueError("The validated FITS product is no longer available on disk.")
            digest = _sha256_file(path)
            if product.sha256 and digest.lower() != str(product.sha256).lower():
                raise ValueError("The cached FITS checksum changed after the product was validated.")

            with fits.open(path, memmap=False) as hdul:
                image_hdu = hdul["IMAGE"] if "IMAGE" in hdul else None
                flags_hdu = hdul["FLAGS"] if "FLAGS" in hdul else None
                variance_hdu = hdul["VARIANCE"] if "VARIANCE" in hdul else None
                if image_hdu is None or flags_hdu is None or variance_hdu is None or "WCS-WAVE" not in hdul:
                    raise ValueError("The product must contain IMAGE, FLAGS, VARIANCE, and WCS-WAVE extensions.")
                image = np.asarray(image_hdu.data, dtype=float)
                flags = np.asarray(flags_hdu.data)
                variance = np.asarray(variance_hdu.data, dtype=float)
                if image.ndim != 2 or flags.shape != image.shape or variance.shape != image.shape:
                    raise ValueError("IMAGE, FLAGS, and VARIANCE must have matching two-dimensional shapes.")
                image_header = image_hdu.header
                image_wcs = WCS(image_header)
                spectral_wcs = WCS(image_header, fobj=hdul, key="W")
                spectral_wcs.sip = None
                coordinate = SkyCoord(ra_deg, dec_deg, unit="deg", frame="icrs")
                x_pixel, y_pixel = image_wcs.world_to_pixel(coordinate)
                if not np.isfinite(x_pixel) or not np.isfinite(y_pixel):
                    raise ValueError("The selected sky coordinate does not map to finite image pixels.")
                # The official SDT photometer uses a 7-to-16 pixel background
                # annulus; require its full support to lie inside the cutout.
                margin = 17.0
                if x_pixel < margin or y_pixel < margin or x_pixel > image.shape[1] - 1 - margin or y_pixel > image.shape[0] - 1 - margin:
                    raise ValueError("The selected position is too close to the cutout edge for the SDT 7-to-16 pixel background annulus.")

                spectral_world = spectral_wcs.pixel_to_world(x_pixel, y_pixel)
                wavelength = float(spectral_world[0].to_value(u.um))
                bandwidth = float(spectral_world[1].to_value(u.um))
                if not np.isfinite(wavelength) or not np.isfinite(bandwidth) or bandwidth <= 0:
                    raise ValueError("WCS-WAVE did not provide finite wavelength and bandwidth values.")
                mjd = float(image_header["MJD-OBS"])
                date_utc = Time(mjd, format="mjd", scale="utc").isot
                detector = int(image_header.get("DETECTOR", 0))
                if detector < 1 or detector > 6:
                    raise ValueError(f"FITS IMAGE detector value {detector} is outside the SPHEREx D1-D6 range.")
                image_unit = u.Unit(image_header.get("BUNIT", ""))
                if image_unit == u.MJy / u.sr:
                    scale_map = sa_mjy_scale_map(
                        band=detector,
                        ny=image.shape[0],
                        nx=image.shape[1],
                        SA_MAPS=sapm_maps,
                        SA_MAPS_interp=sapm_interpolators,
                    )
                    image_for_phot = image * scale_map  # mJy per pixel
                    variance_for_phot = variance * scale_map**2  # mJy squared per pixel
                    flux_unit = "mJy"
                    abmag_scale = 1e-3  # mJy to Jy
                else:
                    image_for_phot = image
                    variance_for_phot = variance
                    flux_unit = str(image_unit)
                    abmag_scale = None

                flag_definitions: dict[str, int] = {}
                for key in flags_hdu.header:
                    try:
                        if key.startswith("MSKN"):
                            flag_definitions[str(flags_hdu.header[key]).upper()] = int(key[4:])
                        elif key.startswith("MP_"):
                            flag_definitions[key[3:].upper()] = int(flags_hdu.header[key])
                    except (TypeError, ValueError):
                        continue
                severe_masks = {
                    name: (flags & (1 << bit)) > 0
                    for name, bit in flag_definitions.items()
                    if name in SDT_BAD_FLAGS
                }
                bad_pixel_mask = np.logical_or.reduce(list(severe_masks.values())) if severe_masks else np.zeros(image.shape, dtype=bool)
                bad_pixel_mask |= ~np.isfinite(image_for_phot) | ~np.isfinite(variance_for_phot) | (variance_for_phot < 0)
                image_masked = np.where(bad_pixel_mask, np.nan, image_for_phot)
                position_xy = [(float(x_pixel), float(y_pixel))]
                radius_list = [aperture_radius_px]
                table = extract_photometry(
                    image_masked,
                    position_xy,
                    radius_list,
                    annulus_rad_in=7,
                    annulus_rad_out=16,
                    mask=bad_pixel_mask,
                )
                flux_key = f"aper_sum_bkgsub_r{aperture_radius_px:.1f}"
                flux_value = float(table[flux_key][0])
                aperture = CircularAperture(position_xy, r=aperture_radius_px)
                variance_result = aperture_photometry(
                    variance_for_phot,
                    aperture,
                    mask=bad_pixel_mask,
                    method="subpixel",
                    subpixels=5,
                )
                variance_sum = float(variance_result["aperture_sum"][0])
                if not np.isfinite(flux_value) or not np.isfinite(variance_sum) or variance_sum < 0:
                    raise ValueError("A finite aperture flux and non-negative aperture variance are required.")
                if abmag_scale is not None:
                    ab_value = float(get_AB(abmag_scale * flux_value))
                    ab_magnitude = ab_value if np.isfinite(ab_value) else None
                else:
                    ab_magnitude = None
                dqf = {
                    name: extract_dqf(flag_mask, position_xy, radius_list)
                    for name, flag_mask in severe_masks.items()
                }
                bad_pixel_count = 0.0
                for name in severe_masks:
                    bad_pixel_count += float(dqf[name][f"aperture_sum_r{aperture_radius_px:.1f}"][0])
                measurement.update({
                    "status": "ok",
                    "date_utc": date_utc,
                    "mjd_utc": mjd,
                    "detector": detector,
                    "product_version": str(hdul[0].header.get("VERSION", "unknown")),
                    "wavelength_um": wavelength,
                    "bandwidth_um": bandwidth,
                    "flux": flux_value,
                    "flux_unit": flux_unit,
                    "flux_error": float(np.sqrt(max(variance_sum, 0.0))),
                    "variance_sum": variance_sum,
                    "ab_magnitude": ab_magnitude,
                    "bad_flagged_aperture_pixels": bad_pixel_count,
                    "known_bad_flags_found": sorted(severe_masks),
                    "x_pixel": float(x_pixel),
                    "y_pixel": float(y_pixel),
                })
        except Exception as exc:
            measurement.update({"status": "failed", "detail": str(exc)})
        output.append(measurement)

    successful = [row for row in output if row.get("status") == "ok"]
    return {
        "status": "ok" if successful and len(successful) == len(output) else "partial" if successful else "failed",
        "source_id": str(source_id),
        "aperture_radius_px": aperture_radius_px,
        "measurements": sorted(successful, key=lambda row: (row["wavelength_um"], row["mjd_utc"])) + [row for row in output if row.get("status") != "ok"],
        "sapm_collection": SDT_SAPM_COLLECTION if tuple(sapm_s3) == SDT_SAPM_S3 else "caller-supplied",
        "sapm_sources": list(sapm_s3),
        "method": "IRSA SPHEREx Source Discovery Tool extract_photometry, extract_dqf, sa_mjy_scale_map, and get_AB components",
        "limitation": "Sparse, time-stamped aperture measurements from the selected products. They are not a variability-corrected continuous spectrum; flux errors omit background-estimation and inter-pixel covariance terms.",
        "reference": "https://caltech-ipac.github.io/irsa-tutorials/spherex-source-discovery-tool-demo/",
    }
