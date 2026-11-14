"""Central, non-user-specific configuration for THE CHAOS THEORY.

Archive endpoints and product conventions were verified against the official
IRSA SPHEREx documentation on 2026-11-13. See ARCHIVE_DOCS below.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent
ASSETS_DIR = ROOT / "assets"
DATA_DIR = ROOT / "data"
REAL_DATA_DIR = DATA_DIR / "real"
DEMO_DATA_DIR = DATA_DIR / "demo"
UPLOAD_DIR = DATA_DIR / "uploads"
CACHE_DIR = DATA_DIR / "cache"
FITS_CACHE_DIR = CACHE_DIR / "fits"
OUTPUT_DIR = ROOT / "outputs"
STATIC_DIR = ROOT / "static"
TESTS_DIR = ROOT / "tests"

for _directory in (
    REAL_DATA_DIR,
    DEMO_DATA_DIR,
    UPLOAD_DIR,
    FITS_CACHE_DIR,
    OUTPUT_DIR,
    STATIC_DIR,
):
    _directory.mkdir(parents=True, exist_ok=True)

APP_TITLE = "THE CHAOS THEORY"
APP_TAGLINE = "Explore the Universe. Detect the Unexpected."
APP_VERSION = "1.1.0-science-integration"

# Verified official service endpoints.
IRSA_TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
IRSA_SIA2 = "https://irsa.ipac.caltech.edu/SIA"
IRSA_CUTOUT_ROOT = "https://irsa.ipac.caltech.edu/"
IRSA_SPHEREX_IBE_ROOT = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
JPL_SBDB_API = "https://ssd-api.jpl.nasa.gov/sbdb.api"
JPL_HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"
SPHEREX_HIPS_ROOT = "https://irsa.ipac.caltech.edu/data/hips/SPHEREx/spherex_allsky_cube"
SPHEREX_HIPS_DATA_ROOT = "https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/hips5_cube/hips_20261026"

# Official SIA2 collections (SPHEREx Data Access, updated 2026-11-12).
SIA2_COLLECTIONS = {
    "DR1 wide": "spherex_dr1",
    "DR1 deep": "spherex_dr1_deep",
    "QR3 wide": "spherex_qr3",
    "QR3 deep": "spherex_qr3_deep",
}

# Documented detector/LVF wavelength ranges and resolving powers.
DETECTOR_BANDS = {
    "SPHEREx-D1": {"range_um": (0.75, 1.09), "R": 39, "label": "D1 · 0.75–1.09 µm"},
    "SPHEREx-D2": {"range_um": (1.10, 1.62), "R": 41, "label": "D2 · 1.10–1.62 µm"},
    "SPHEREx-D3": {"range_um": (1.63, 2.41), "R": 41, "label": "D3 · 1.63–2.41 µm"},
    "SPHEREx-D4": {"range_um": (2.42, 3.82), "R": 35, "label": "D4 · 2.42–3.82 µm"},
    "SPHEREx-D5": {"range_um": (3.83, 4.41), "R": 112, "label": "D5 · 3.83–4.41 µm"},
    "SPHEREx-D6": {"range_um": (4.42, 5.00), "R": 128, "label": "D6 · 4.42–5.00 µm"},
}

# SPHEREx FLAGS bit names read dynamically from the FITS header.
# SOURCE and FULLSAMPLE are informational, not defects; masking them would
# remove almost every pixel in a real cutout.
BAD_FLAG_NAMES = {
    "TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "DICHROIC",
    "MISSING_DATA", "HOT", "COLD", "PHANMISS", "NONLINEAR", "PERSIST",
    "OUTLIER", "CROSSTALK", "GHOST", "GHOST_FPA", "GHOST_EXT", "STREAK",
    "BLOOM", "SNOWBALL", "HALO",
}
INFO_FLAG_NAMES = {"SOURCE", "FULLSAMPLE"}

DATA_MODES = {
    "EMPTY": "NO DATA LOADED",
    "LIVE": "LIVE IRSA SPHEREx",
    "CACHED_REAL": "REAL SPHEREx · CACHED IRSA CUTOUTS",
    "LOCAL": "LOCAL FITS UPLOAD · UNVERIFIED PROVENANCE",
    "DEMO": "DEMONSTRATION DATA · SYNTHETIC",
}

# Verified M101 field from the official IRSA cutout tutorial. The cached pair
# is documented in data/real/manifest.json and was retrieved from live IRSA.
DEFAULT_QUERY = {
    "ra_deg": 210.80227,
    "dec_deg": 54.34895,
    "radius_deg": 0.05,
    "size_deg": 0.10,
    "start_date": "2026-05-01",
    "end_date": "2026-06-30",
    "band": "SPHEREx-D2",
    "release": "all",
    "survey": "wide",
    "max_rows": 120,
}

PROCESSING_DEFAULTS = {
    "target_reference": "epoch1",
    "wavelength_tolerance_bandwidth_fraction": 0.25,
    "min_overlap_fraction": 0.55,
    "background_sigma": 3.0,
    "background_iterations": 5,
    "fit_background_plane": True,
    "normalization_mode": "background + robust linear scaling",
    "robust_fit_samples": 5000,
    "detection_sigma": 5.0,
    "min_distance_px": 3.0,
    "border_px": 3,
    "dipole_distance_px": 4.0,
    "aperture_radius_px": 2.0,
    "annulus_inner_px": 4.0,
    "annulus_outer_px": 6.0,
    "max_candidates": 40,
    "max_flag_fraction": 0.20,
    "crossmatch_enabled": True,
    "gaia_crossmatch_enabled": True,
    "crossmatch_catalogue": "2MASS Point Source Catalog (IRSA TAP fp_psc)",
    "crossmatch_radius_arcsec": 6.0,
    "catalogue_cone_radius_deg": 0.08,
    "variance_sigma_clip_fraction": 0.2,
    "variance_sigma_clip_max": 5.0,
}

# Memory and network safety limits for limited hackathon hardware.
MAX_UPLOAD_BYTES = 200 * 1024 * 1024
MAX_DOWNLOAD_BYTES = 120 * 1024 * 1024
MAX_FITS_PIXELS = 4_000_000
DEFAULT_CROP_SIDE = 512
REQUEST_TIMEOUT_SECONDS = 120
ARCHIVE_TIMEOUT_SECONDS = 90
CACHE_TTL_SECONDS = 1800

ARCHIVE_DOCS = [
    {
        "title": "SPHEREx Data Access — IRSA archive documentation",
        "url": "https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/",
    },
    {
        "title": "SPHEREx Data Products — MEF extension and FLAGS definitions",
        "url": "https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/",
    },
    {
        "title": "IRSA tutorial — SPHEREx spectral image cutouts via TAP",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-cutouts/",
    },
    {
        "title": "IRSA tutorial — Introduction to SPHEREx spectral images",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/",
    },
    {
        "title": "IRSA TAP documentation",
        "url": "https://irsa.ipac.caltech.edu/docs/program_interface/TAP.html",
    },
    {
        "title": "IRSA cutout service parameters",
        "url": "https://irsa.ipac.caltech.edu/ibe/cutouts.html",
    },
    {
        "title": "NASA SPHEREx mission",
        "url": "https://science.nasa.gov/mission/spherex/",
    },
    {
        "title": "ESA Gaia DR3 gaia_source data model",
        "url": "https://gea.esac.esa.int/archive/documentation/GDR3/Gaia_archive/chap_datamodel/sec_dm_main_source_catalogue/ssec_dm_gaia_source.html",
    },
    {
        "title": "ESA Gaia DR3 astrometric propagation and error propagation",
        "url": "https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html",
    },
    {
        "title": "IRSA SPHEREx ePSF tutorial",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/",
    },
    {
        "title": "IRSA SPHEREx Source Discovery Tool tutorial",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-source-discovery-tool-demo/",
    },
    {
        "title": "JPL Small-Body Database API",
        "url": "https://ssd-api.jpl.nasa.gov/doc/sbdb.html",
    },
    {
        "title": "JPL Horizons API",
        "url": "https://ssd-api.jpl.nasa.gov/doc/horizons.html",
    },
    {
        "title": "IRSA SPHEREx Mosaic Cube tutorial",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-mosaic/",
    },
    {
        "title": "IRSA SPHEREx HiPS All-Sky Cube tutorial",
        "url": "https://caltech-ipac.github.io/irsa-tutorials/spherex-hipscube/",
    },
]

NAV_ITEMS = [
    ("Overview", "pages/01_Mission_Control.py", "01", "\u25ce"),
    ("Explore Sky", "pages/02_Sky_Explorer.py", "02", "\u2316"),
    ("Observations", "pages/03_Observation_Explorer.py", "03", "\u2637"),
    ("Image Comparison", "pages/04_Cosmic_Image_Lab.py", "04", "\u25a3"),
    ("Detection Engine", "pages/05_Detection_Engine.py", "05", "\u2699"),
    ("Candidate Analysis", "pages/06_Candidate_Investigation.py", "06", "\u2727"),
    ("Data Explorer", "pages/08_Data_Explorer.py", "08", "\u2261"),
    ("Methodology", "pages/07_Methodology.py", "07", "\u2193"),
    ("About", "pages/09_About.py", "09", "\u24d8"),
]
