"""Optional adapter for the official IRSA SPHEREx Source Discovery Tool.

The core app remains usable without the package or a SExtractor executable.
This adapter uses the IRSA package's packaged SExtractor configuration files
and invokes the installed executable directly, keeping this integration
separate from the Photutils residual detector.
"""
from __future__ import annotations

import importlib.util
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.table import Table

from modules.models import FitsProduct


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
        "interpretation": "Spatial agreement diagnostic only; the two extraction paths use different algorithms and neither count is a calibrated completeness estimate.",
    }
