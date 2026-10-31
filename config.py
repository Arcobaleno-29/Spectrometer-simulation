import numpy as np
import math
from pathlib import Path
from dataclasses import dataclass

@dataclass
class OpticsConfig:
    slit_width_um: float 
    slit_length_mm: float
    collimator_focal_length_mm: float
    collimator_distance_mm: float
    focusing_focal_length_mm: float
    lens_separation_mm: float
    detector_distance_mm: float
    pupil_diameter_mm: float
    poppy_npix: int

    def __post_init__(self):

        # Positive values check
        positive_parameters = {
            "slit_width_um": self.slit_width_um,
            "slit_length_mm": self.slit_length_mm,
            "collimator_distance_mm": self.collimator_distance_mm,
            "lens_separation_mm": self.lens_separation_mm,
            "detector_distance_mm": self.detector_distance_mm,
            "pupil_diameter_mm": self.pupil_diameter_mm,
        }

        for name, value in positive_parameters.items():
            if value <= 0:
                raise ValueError(f"{name} must be positive.")

        # Sampling needs an integer check 
        if (
            isinstance(self.poppy_npix, bool)
            or not isinstance(self.poppy_npix, int)
            or self.poppy_npix <= 0
        ):
            raise ValueError("poppy_npix must be a positive integer.")

@dataclass
class RGBComponent:
    name: str
    center_nm: float
    fwhm_nm: float
    amplitude: float = 1.0

    def __post_init__(self):
        if not self.name:
            raise ValueError("name must not be empty")
        if not math.isfinite(self.center_nm) or self.center_nm <= 0:
            raise ValueError("center_nm must be finite and positive")
        if not math.isfinite(self.fwhm_nm) or self.fwhm_nm <= 0:
            raise ValueError("fwhm_nm must be finite and positive")
        if not math.isfinite(self.amplitude) or self.amplitude < 0:
            raise ValueError("amplitude must be finite and non-negative")
        
@dataclass
class DetectorConfig:
    width_px: int
    height_px: int
    pixel_size_um: float
    psf_fov_pixels: int
    bayer_pattern: str = "RGGB"

    def __post_init__(self):
        integer_parameters = {
            "width_px": self.width_px,
            "height_px": self.height_px,
            "psf_fov_pixels": self.psf_fov_pixels
        }

        for name, value in integer_parameters.items():
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or value <= 0
            ):
                raise ValueError(
                    f"{name} must be a positive integer."
                )

        if (
            not math.isfinite(self.pixel_size_um)
            or self.pixel_size_um <= 0
        ):
            raise ValueError(
                "pixel_size_um must be finite and positive."
            )

        valid_patterns = {"RGGB", "BGGR", "GRBG", "GBRG"}

        if self.bayer_pattern not in valid_patterns:
            raise ValueError(
                f"bayer_pattern must be one of {sorted(valid_patterns)}."
            )
@dataclass
class GratingConfig:
    lines_per_mm: float
    diffraction_order: int
    incidence_angle_deg: float

    def __post_init__(self):
        if not math.isfinite(self.lines_per_mm) or self.lines_per_mm <= 0:
            raise ValueError("lines_per_mm must be finite and positive")

        if (
            isinstance(self.diffraction_order, bool)
            or not isinstance(self.diffraction_order, int)
            or self.diffraction_order < 1
        ):
            raise ValueError("diffraction_order must be a positive integer")

        if not math.isfinite(self.incidence_angle_deg):
            raise ValueError("incidence_angle_deg must be finite")

@dataclass
class SpectrumConfig:
    min_wavelength_nm: float
    max_wavelength_nm: float
    samples: int
    center_wavelength_nm: float

    def __post_init__(self):
        positive_parameters = {
            "min_wavelength_nm": self.min_wavelength_nm,
            "max_wavelength_nm": self.max_wavelength_nm,
            "center_wavelength_nm": self.center_wavelength_nm,
        }

        for name, value in positive_parameters.items():
            if value <= 0:
                raise ValueError(f"{name} must be positive.")

        integer_parameters = {"samples": self.samples}

        for name, value in integer_parameters.items():
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or value <= 0
            ):
                raise ValueError(
                    f"{name} must be a positive integer."
                )
@dataclass
class AnalysisConfig:
    resolution_test_wavelength_nm: float = 550.0
    resolution_curve_points: int = 21
    point_source_radius_um: float = 5.0

@dataclass
class SpectrographConfig:
    optics: OpticsConfig
    grating: GratingConfig
    detector: DetectorConfig
    spectrum: SpectrumConfig
    rgb_source: tuple[RGBComponent, ...]
    analysis: AnalysisConfig
    qe_data_path: Path


def validate_config(cfg: SpectrographConfig, qe) -> None:
    """Check that the instrument and measured QE support the simulated band."""
    spectrum = cfg.spectrum
    if spectrum.min_wavelength_nm >= spectrum.max_wavelength_nm:
        raise ValueError("spectrum minimum must be below its maximum")
    if not (
        spectrum.min_wavelength_nm
        <= spectrum.center_wavelength_nm
        <= spectrum.max_wavelength_nm
    ):
        raise ValueError("spectrum center must be within the wavelength range")

    qe_wavelengths = np.asarray(qe.wavelengths_nm, dtype=float)
    if (
        qe_wavelengths.ndim != 1
        or qe_wavelengths.size < 2
        or not np.all(np.isfinite(qe_wavelengths))
        or not np.all(np.diff(qe_wavelengths) > 0)
    ):
        raise ValueError("QE wavelengths must be finite and strictly increasing")
    if (
        spectrum.min_wavelength_nm < qe_wavelengths[0]
        or spectrum.max_wavelength_nm > qe_wavelengths[-1]
    ):
        raise ValueError("Simulation wavelengths exceed the measured QE range")
    for channel in (qe.blue, qe.green, qe.red):
        values = np.asarray(channel, dtype=float)
        if values.shape != qe_wavelengths.shape:
            raise ValueError("Every QE channel must match the wavelength table")
        if not np.all(np.isfinite(values)) or np.any((values < 0) | (values > 1)):
            raise ValueError("QE values must be finite fractions between 0 and 1")

    spacing_nm = 1e6 / cfg.grating.lines_per_mm
    alpha = math.radians(cfg.grating.incidence_angle_deg)
    for wavelength_nm in (spectrum.min_wavelength_nm, spectrum.max_wavelength_nm):
        argument = (
            cfg.grating.diffraction_order * wavelength_nm / spacing_nm
            - math.sin(alpha)
        )
        if abs(argument) > 1:
            raise ValueError(
                "The selected grating/order cannot diffract the full wavelength range"
            )


def default_config() -> SpectrographConfig:
    return SpectrographConfig (
        optics=OpticsConfig(
            slit_width_um=20.0,
            slit_length_mm=3.0,
            collimator_focal_length_mm=75.0,
            collimator_distance_mm=75.0,
            focusing_focal_length_mm=30.0,
            lens_separation_mm=30.0,
            detector_distance_mm=30.0,
            pupil_diameter_mm=25.4,
            poppy_npix=512,
            ),
        grating=GratingConfig(
            lines_per_mm=300.0, 
            diffraction_order=1, 
            incidence_angle_deg=0.0
            ),
        detector=DetectorConfig(
            width_px=1440, 
            height_px=1080, 
            pixel_size_um=3.45,
            psf_fov_pixels= 512, 
            bayer_pattern="RGGB"),
        spectrum=SpectrumConfig(
            min_wavelength_nm=400.0, 
            max_wavelength_nm=700.0, 
            samples=601, 
            center_wavelength_nm=551.7),
        rgb_source=(
        RGBComponent("blue", center_nm=450.0, fwhm_nm=10.0),
        RGBComponent("green", center_nm=530.0, fwhm_nm=10.0),
        RGBComponent("red", center_nm=620.0, fwhm_nm=10.0),
        ),
        analysis=AnalysisConfig(),
        qe_data_path=Path(__file__).resolve().parent
        / "data/cs165cu_quantum_efficiency.csv",
    )
