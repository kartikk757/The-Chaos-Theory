"""FITS acquisition, local upload handling and validated product loading."""
from __future__ import annotations

import hashlib
import json
import re
import shutil
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import requests
import streamlit as st
from astropy.io import fits
from astropy.time import Time

import config
from modules import coordinates, fits_processor, provenance
from modules.network_guard import wait_for_service
from modules.models import FitsProduct


class DataLoadError(RuntimeError):
    pass


def _resolve_path(path: str | Path) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = config.ROOT / candidate
    return candidate


@st.cache_data(ttl=3600, max_entries=32, show_spinner=False)
def fetch_cutout_path(url: str, force: bool = False, max_bytes: int = config.MAX_DOWNLOAD_BYTES) -> dict:
    """Download and cache a documented IRSA cutout URL. Never fabricate a file."""
    if not url.startswith("https://irsa.ipac.caltech.edu/"):
        raise DataLoadError("Only verified NASA/IPAC IRSA URLs are accepted for remote FITS retrieval")
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:24]
    path = config.FITS_CACHE_DIR / f"{digest}.fits"
    sidecar = config.FITS_CACHE_DIR / f"{digest}.json"
    if path.exists() and not force and path.stat().st_size > 1000:
        with path.open("rb") as handle:
            if handle.read(8).startswith(b"SIMPLE  "):
                meta = json.loads(sidecar.read_text(encoding="utf-8")) if sidecar.exists() else {}
                return {"path": str(path), "sha256": meta.get("sha256", provenance.sha256_file(path)), "bytes": path.stat().st_size, "url": url, "cache": "local"}
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            wait_for_service("irsa_cutout", 0.5)
            response = requests.get(url, stream=True, timeout=config.REQUEST_TIMEOUT_SECONDS)
            response.raise_for_status()
            declared = int(response.headers.get("Content-Length") or 0)
            if declared > max_bytes:
                raise DataLoadError(f"Remote cutout declares {declared:,} bytes, above the {max_bytes:,} byte safety limit")
            temp = path.with_suffix(".partial")
            total = 0
            sha = hashlib.sha256()
            with temp.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=256 * 1024):
                    if not chunk:
                        continue
                    total += len(chunk)
                    if total > max_bytes:
                        raise DataLoadError(f"Remote cutout exceeded the {max_bytes:,} byte safety limit")
                    handle.write(chunk)
                    sha.update(chunk)
            if total < 1000 or not temp.open("rb").read(8).startswith(b"SIMPLE  "):
                temp.unlink(missing_ok=True)
                raise DataLoadError("IRSA response was not a valid FITS file (missing SIMPLE header)")
            temp.replace(path)
            checksum = sha.hexdigest()
            meta = {
                "url": url,
                "sha256": checksum,
                "bytes": total,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "content_type": response.headers.get("Content-Type", ""),
            }
            sidecar.write_text(json.dumps(meta, indent=2), encoding="utf-8")
            return {"path": str(path), "sha256": checksum, "bytes": total, "url": url, "cache": "network", "retrieved_utc": meta["retrieved_utc"]}
        except Exception as exc:
            last_error = exc
            time.sleep(1.0 + attempt)
    raise DataLoadError(f"IRSA cutout retrieval failed after 3 attempts: {last_error}")


def _load_product_impl(path: str, file_sha256: str, center_ra: float | None, center_dec: float | None, record_json: str) -> FitsProduct:
    record = json.loads(record_json)
    resolved = _resolve_path(path)
    warnings_out: list[str] = []
    with fits.open(resolved, memmap=False) as hdul:
        primary_header = fits_processor.header_dict(hdul[0].header)
        extension_names = [str(hdu.name) for hdu in hdul]
        image_hdu = fits_processor.find_image_hdu(hdul)
        image_header = image_hdu.header
        image = np.asarray(image_hdu.data, dtype=float)
        flags_hdu = fits_processor.find_extension(hdul, "FLAGS")
        variance_hdu = fits_processor.find_extension(hdul, "VARIANCE")
        zodi_hdu = fits_processor.find_extension(hdul, "ZODI")
        flags = np.asarray(flags_hdu.data) if flags_hdu is not None and getattr(flags_hdu, "data", None) is not None else None
        variance = np.asarray(variance_hdu.data, dtype=float) if variance_hdu is not None and getattr(variance_hdu, "data", None) is not None else None
        zodi = np.asarray(zodi_hdu.data, dtype=float) if zodi_hdu is not None and getattr(zodi_hdu, "data", None) is not None else None
        header = fits_processor.header_dict(image_header)
        bunit = str(image_header.get("BUNIT", "") or "")
        definitions = fits_processor.parse_flag_definitions(flags_hdu.header if flags_hdu is not None else None)
        bad_mask, source_mask, flag_diag = fits_processor.build_flag_masks(flags, definitions)
        wcs, wcs_warnings = fits_processor.build_spatial_wcs(image_header)
        warnings_out.extend(wcs_warnings)
        # The delivered cutout is centred on the requested sky position, which is NOT
        # necessarily the observation footprint centre recorded in the archive. The
        # LVF central wavelength varies across the focal plane, so it must be sampled
        # at a pixel that exists in this file: the geometric centre of the array.
        image_center_px = ((image.shape[1] - 1) / 2.0, (image.shape[0] - 1) / 2.0)
        center_wl, center_bw, spectral_warnings = fits_processor.spectral_measurement(
            hdul, image_header, wcs, image_center_px[0], image_center_px[1]
        )
        warnings_out.extend(spectral_warnings)
        center_world = fits_processor.center_sky(wcs, image.shape)
        coordinate_info = coordinates.coordinate_context(record)
        requested_sky = coordinate_info.get("requested_field_center_icrs_deg")
        requested_px = None
        if wcs is not None and requested_sky is not None:
            try:
                from astropy.coordinates import SkyCoord as _SkyCoord

                x_req, y_req = wcs.world_to_pixel(
                    _SkyCoord(float(requested_sky[0]), float(requested_sky[1]), unit="deg", frame="icrs")
                )
                if -0.5 <= x_req <= image.shape[1] - 0.5 and -0.5 <= y_req <= image.shape[0] - 0.5:
                    requested_px = (float(x_req), float(y_req))
                else:
                    warnings_out.append(
                        "The requested field centre lies outside the delivered cutout; it is retained as query provenance, "
                        "while cutout geometry is measured from the FITS WCS."
                    )
            except Exception:
                requested_px = None
        center_px = requested_px or image_center_px
        scale = fits_processor.pixel_scale_arcsec(wcs)
        arrays = {"image": image, "flags": flags, "variance": variance, "zodi": zodi, "bad_mask": bad_mask, "source_mask": source_mask}
        cropped, cropped_wcs, cropped_center, crop_warnings = fits_processor.crop_to_limit(arrays, wcs)
        if crop_warnings:
            warnings_out.extend(crop_warnings)
            # The WCS-WAVE lookup table is tied to original detector coordinates; do not
            # present an unadjusted wavelength for a cropped local file.
            center_wl, center_bw = None, None
            warnings_out.append("Central wavelength was not recomputed after cropping; spectral compatibility is unverified.")
            wcs = cropped_wcs
            center_px = cropped_center
        original_shape = image.shape
        image = cropped["image"]
        flags = cropped["flags"]
        variance = cropped["variance"]
        zodi = cropped["zodi"]
        bad_mask = np.asarray(cropped["bad_mask"], dtype=bool)
        source_mask = np.asarray(cropped["source_mask"], dtype=bool)
    if flags is not None and np.asarray(flags).shape != image.shape:
        flags = None
        warnings_out.append("FLAGS shape became incompatible after cropping and was discarded.")
    if variance is not None and np.asarray(variance).shape != image.shape:
        variance = None
        warnings_out.append("VARIANCE shape became incompatible after cropping and was discarded.")
    if zodi is not None and np.asarray(zodi).shape != image.shape:
        zodi = None
        warnings_out.append("ZODI shape became incompatible after cropping and was discarded.")
    if bad_mask.shape != image.shape:
        bad_mask = np.zeros(image.shape, dtype=bool)
    if source_mask.shape != image.shape:
        source_mask = np.zeros(image.shape, dtype=bool)
    if flags is None:
        definitions = {}
        flag_diag = {"available": False}
    if not np.isfinite(image).any():
        raise DataLoadError("FITS IMAGE contains no finite pixels")
    if np.mean(np.isfinite(image)) < 0.5:
        warnings_out.append("More than half of IMAGE pixels are non-finite; results will be heavily masked.")
    if flag_diag.get("available") and flag_diag.get("bad_pixel_fraction", 0) > 0.5:
        warnings_out.append("More than half of pixels carry severe processing flags; check the FLAGS definitions before interpreting residuals.")
    validation = fits_processor.validate_arrays(
        image, flags, variance, zodi, wcs, bunit, bad_mask, definitions, center_wl, scale
    )
    metadata = {
        "observation": record,
        "file_path": str(resolved),
        "file_sha256": file_sha256,
        "file_size_bytes": resolved.stat().st_size,
        "access_url": record.get("cutout_url") or record.get("full_fits_url") or "",
        "flag_diagnostics": flag_diag,
        "wavelength_sample_pixel": [image_center_px[0], image_center_px[1]],
        "wavelength_sample_note": "WCS-WAVE evaluated at the geometric centre of the delivered cutout array.",
        "requested_position_pixel": list(requested_px) if requested_px else None,
        "record_position_icrs_deg": coordinate_info.get("record_position_icrs_deg"),
        "record_position_role": coordinate_info.get("record_position_role"),
        "record_position_source": coordinate_info.get("record_position_source"),
        "requested_field_center_icrs_deg": requested_sky,
        "requested_field_center_source": coordinate_info.get("requested_field_center_source"),
        "coordinate_frame": coordinate_info.get("coordinate_frame", "ICRS"),
        "cutout_center_sky": list(center_world) if center_world else None,
        "crop_origin_pixel": [
            float((original_shape[1] - image.shape[1]) / 2.0) if crop_warnings else 0.0,
            float((original_shape[0] - image.shape[0]) / 2.0) if crop_warnings else 0.0,
        ],
        "extension_names": extension_names,
        "primary_version": primary_header.get("VERSION", ""),
        "loaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return FitsProduct(
        product_id=str(record.get("id") or resolved.stem),
        path=str(resolved),
        origin=str(record.get("origin") or "FITS file"),
        data_mode=str(record.get("data_mode") or config.DATA_MODES["LOCAL"]),
        image=image,
        flags=flags,
        variance=variance,
        zodi=zodi,
        wcs=wcs,
        header=header,
        primary_header=primary_header,
        extension_names=extension_names,
        flag_definitions=definitions,
        bad_flag_mask=bad_mask,
        source_flag_mask=source_mask,
        bunit=bunit,
        shape=tuple(int(v) for v in image.shape),
        center_wavelength_um=center_wl,
        center_bandwidth_um=center_bw,
        center_pixel=(float(center_px[0]), float(center_px[1])),
        center_sky=(float(center_world[0]), float(center_world[1])) if center_world else None,
        pixel_scale_arcsec=scale,
        file_size_bytes=int(resolved.stat().st_size),
        sha256=file_sha256,
        metadata=metadata,
        warnings=warnings_out + list(record.get("warnings", [])),
        validation=validation,
        is_real_spherex=bool(record.get("is_real_spherex", False)),
        is_demo=bool(record.get("is_demo", False)),
    )


@st.cache_data(max_entries=16, show_spinner=False)
def _load_fits_product_cached(path: str, mtime: float, file_sha256: str, center_ra: float | None, center_dec: float | None, record_json: str) -> FitsProduct:
    return _load_product_impl(path, file_sha256, center_ra, center_dec, record_json)


def load_fits_product(path: str | Path, record: dict, center_ra: float | None = None, center_dec: float | None = None) -> FitsProduct:
    resolved = _resolve_path(path)
    if not resolved.exists():
        raise DataLoadError(f"FITS file does not exist: {resolved}")
    if resolved.stat().st_size > config.MAX_DOWNLOAD_BYTES:
        raise DataLoadError(f"FITS file is larger than the {config.MAX_DOWNLOAD_BYTES:,} byte safety limit")
    checksum = provenance.sha256_file(resolved)
    record_json = json.dumps(provenance.safe_json(record), sort_keys=True)
    return _load_fits_product_cached(
        str(resolved), resolved.stat().st_mtime, checksum, center_ra, center_dec, record_json
    )


def load_observation_product(record: dict, force_download: bool = False) -> FitsProduct:
    cached_path = record.get("cached_path")
    if cached_path:
        path = _resolve_path(cached_path)
        if not path.exists():
            raise DataLoadError(
                f"Cached FITS file does not exist: {path}. Re-retrieve this cutout from IRSA or choose another "
                "observation; the application will not substitute a different file."
            )
    elif record.get("cutout_url"):
        fetched = fetch_cutout_path(record["cutout_url"], force=force_download)
        path = Path(fetched["path"])
    else:
        raise DataLoadError("This record has neither a cached FITS path nor a documented IRSA cutout URL")
    return load_fits_product(path, record, center_ra=record.get("ra_deg"), center_dec=record.get("dec_deg"))


def save_uploaded_fits(uploaded_file: Any) -> dict:
    name = str(getattr(uploaded_file, "name", "upload.fits"))
    suffix = Path(name).suffix.lower()
    if suffix not in {".fits", ".fit", ".fts", ".gz"} and not name.lower().endswith(".fits.gz"):
        raise DataLoadError("Only .fits, .fit, .fts or .fits.gz files are accepted")
    payload = bytes(uploaded_file.getbuffer())
    if len(payload) > config.MAX_UPLOAD_BYTES:
        raise DataLoadError(f"Uploaded FITS exceeds the {config.MAX_UPLOAD_BYTES:,} byte safety limit")
    if not payload.startswith(b"SIMPLE  ") and suffix != ".gz":
        raise DataLoadError("Uploaded file does not begin with a FITS SIMPLE header")
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", Path(name).name)[:100]
    path = config.UPLOAD_DIR / f"{uuid.uuid4().hex[:8]}_{safe}"
    path.write_bytes(payload)
    return {
        "path": str(path),
        "filename": name,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "uploaded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def make_local_observation(upload: dict, product: FitsProduct, index: int) -> dict:
    header = product.header
    date_obs = str(header.get("DATE-OBS") or header.get("DATEOBS") or header.get("DATE") or "").strip()
    mjd = None
    date_utc = None
    if date_obs:
        try:
            mjd = float(Time(date_obs, format="isot").mjd)
            date_utc = Time(mjd, format="mjd").isot[:19]
        except Exception:
            pass
    band = str(header.get("BAND") or header.get("FILTER") or header.get("DETECTOR") or "").strip()
    if header.get("DETECTOR") and not band.startswith("SPHEREx"):
        band = f"local detector {band}"
    if not band:
        band = "unknown (local FITS)"
    quality = str(header.get("QUALITY") or header.get("QUALFLAG") or "").strip()
    stable = hashlib.sha256((upload["path"] + upload["sha256"]).encode()).hexdigest()[:12].upper()
    center = product.center_sky or (None, None)
    warnings = list(product.warnings)
    if not date_obs:
        warnings.append("No DATE-OBS/DATE observation timestamp was found in the FITS header.")
    if product.wcs is None:
        warnings.append("No celestial WCS; direct comparison will be rejected.")
    return {
        "id": f"LOCAL-{stable}",
        "label": f"Local FITS {index} · {Path(upload['filename']).name}",
        "origin": "User-supplied local FITS upload (provenance unverified)",
        "data_mode": config.DATA_MODES["LOCAL"],
        "ra_deg": center[0],
        "dec_deg": center[1],
        "date_utc": date_utc,
        "mjd": mjd,
        "band": band,
        "detector": str(header.get("DETECTOR", "") or ""),
        "release": "local",
        "obs_id": str(header.get("OBSID") or header.get("OBS_ID") or ""),
        "cutout_url": "",
        "full_fits_url": "",
        "cached_path": upload["path"],
        "content_length": upload["bytes"],
        "calibration_level": None,
        "quality_flag": quality,
        "provenance_version": str(product.primary_header.get("VERSION", "") or ""),
        "exposure_s": float(header["EXPTIME"]) if str(header.get("EXPTIME", "")).replace(".", "", 1).isdigit() else None,
        "energy_bounds_lower_m": None,
        "energy_bounds_upper_m": None,
        "product_type": "local FITS",
        "content_type": "image/fits",
        "plane_id": "",
        "is_real_spherex": False,
        "is_demo": False,
        "central_wavelength_um": product.center_wavelength_um,
        "central_bandwidth_um": product.center_bandwidth_um,
        "coordinates_source": "FITS WCS image centre" if product.wcs is not None else "not available",
        "coordinates_role": "fits_wcs_image_center" if product.wcs is not None else "unavailable",
        "metadata": {"upload": upload, "header_keys": sorted(header.keys()), "sha256": upload["sha256"]},
        "warnings": warnings,
    }


def load_cached_real_bundle() -> tuple[list[dict], dict]:
    manifest = provenance.read_json(config.REAL_DATA_DIR / "manifest.json")
    records: list[dict] = []
    for entry in manifest.get("observations", []):
        record = dict(entry)
        path = _resolve_path(record.get("cached_path", ""))
        if not path.exists():
            # Source repositories can omit the multi-megabyte FITS payloads
            # while retaining their exact IRSA URLs and expected hashes. Fetch
            # only that named product, validate it, and cache it at the
            # manifest path; never replace it with a different observation.
            download_record = dict(record)
            download_record.pop("cached_path", None)
            if not download_record.get("cutout_url") and not download_record.get("full_fits_url"):
                raise DataLoadError(f"Cached real cutout is missing and has no IRSA URL: {path}")
            try:
                product = load_observation_product(download_record)
            except Exception as exc:
                raise DataLoadError(f"Cached real cutout is missing and its IRSA retrieval failed: {exc}") from exc
            expected = str(record.get("derived_sha256") or "").lower()
            source_path = Path(product.path)
            source_sha256 = provenance.sha256_file(source_path)
            if not expected or source_sha256 != expected:
                raise DataLoadError(
                    f"IRSA returned bytes that do not match the manifest checksum for {path.name}; "
                    "the cached reference product was not replaced."
                )
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary_path = path.with_name(path.name + ".partial")
            try:
                shutil.copyfile(source_path, temporary_path)
                if provenance.sha256_file(temporary_path) != expected:
                    raise DataLoadError(f"Copied real cutout checksum failed verification: {temporary_path.name}")
                temporary_path.replace(path)
            finally:
                temporary_path.unlink(missing_ok=True)
        checksum = provenance.sha256_file(path)
        if checksum != record.get("derived_sha256"):
            raise DataLoadError(
                f"Cached real cutout checksum differs from the manifest: {path}. "
                "Remove or re-download this file before scientific use."
            )
        record["cached_path"] = str(path)
        source = str(record.get("coordinates_source") or "")
        record.setdefault(
            "coordinates_role",
            "parent_observation_footprint_center" if "SIA2" in source or "footprint" in source.lower()
            else "requested_field_center" if "requested" in source.lower()
            else "recorded_observation_position",
        )
        # Older manifests store the requested field only in the IRSA cutout URL.
        # Preserve it explicitly for future offline loads and provenance exports.
        coordinate_info = coordinates.coordinate_context(record)
        requested = coordinate_info.get("requested_field_center_icrs_deg")
        if requested:
            coordinates.attach_requested_field_center(record, requested[0], requested[1])
        record["data_mode"] = config.DATA_MODES["CACHED_REAL"]
        record["is_real_spherex"] = True
        record["is_demo"] = False
        records.append(record)
    return records, manifest
