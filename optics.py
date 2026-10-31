import numpy as np
import astropy.units as u
import poppy
import poppy.fresnel as fresnel
from typing import Literal

from config import SpectrographConfig


MappingConvention = Literal["angle_difference", "tangent_difference"]


def _diffraction_angle_nm(wavelength_nm, config: SpectrographConfig):
    """Grating diffraction angle in radians for numeric wavelengths in nm."""
    spacing_nm = 1e6 / config.grating.lines_per_mm
    alpha = np.deg2rad(config.grating.incidence_angle_deg)
    argument = (
        config.grating.diffraction_order
        * np.asarray(wavelength_nm, dtype=float)
        / spacing_nm
        - np.sin(alpha)
    )
    if np.any(~np.isfinite(argument)) or np.any(np.abs(argument) > 1):
        raise ValueError("No real diffraction angle for this wavelength")
    return np.arcsin(argument)


def diffraction_angle(wavelength: u.Quantity, config: SpectrographConfig):
    """Return the grating diffraction angle in radians."""
    return _diffraction_angle_nm(wavelength.to_value(u.nm), config)


def wavelength_to_detector_offset_mm(
    wavelength_nm,
    config: SpectrographConfig,
    *,
    convention: MappingConvention,
):
    """Map wavelength in nm to position in mm relative to the center ray.

    ``angle_difference`` uses x = f tan(beta - beta_center), as in the
    spectrograph simulation. ``tangent_difference`` uses
    x = f (tan(beta) - tan(beta_center)), as in the optimizer notebook.
    Both accept scalar or array wavelengths and put the center ray at x=0.
    """
    beta = _diffraction_angle_nm(wavelength_nm, config)
    beta_center = _diffraction_angle_nm(
        config.spectrum.center_wavelength_nm, config
    )
    focus_mm = config.optics.focusing_focal_length_mm
    if convention == "angle_difference":
        return focus_mm * np.tan(beta - beta_center)
    if convention == "tangent_difference":
        return focus_mm * (np.tan(beta) - np.tan(beta_center))
    raise ValueError(f"Unknown detector mapping convention: {convention}")


def detector_offset_mm_to_wavelength_nm(
    offset_mm,
    config: SpectrographConfig,
    *,
    convention: MappingConvention,
):
    """Invert a detector offset in mm using the specified mapping convention."""
    offset_mm = np.asarray(offset_mm, dtype=float)
    if np.any(~np.isfinite(offset_mm)):
        raise ValueError("Detector offsets must be finite")
    beta_center = _diffraction_angle_nm(
        config.spectrum.center_wavelength_nm, config
    )
    scaled_offset = offset_mm / config.optics.focusing_focal_length_mm
    if convention == "angle_difference":
        beta = beta_center + np.arctan(scaled_offset)
    elif convention == "tangent_difference":
        beta = np.arctan(np.tan(beta_center) + scaled_offset)
    else:
        raise ValueError(f"Unknown detector mapping convention: {convention}")

    spacing_nm = 1e6 / config.grating.lines_per_mm
    alpha = np.deg2rad(config.grating.incidence_angle_deg)
    return (
        spacing_nm
        * (np.sin(alpha) + np.sin(beta))
        / config.grating.diffraction_order
    )


def detector_column_to_offset_mm(column, config: SpectrographConfig):
    """Convert a detector column to mm relative to the array center."""
    return (
        (np.asarray(column, dtype=float) - config.detector.width_px / 2.0)
        * config.detector.pixel_size_um
        / 1000.0
    )


def detector_offset_mm_to_column(offset_mm, config: SpectrographConfig):
    """Convert an offset in mm to a fractional detector column."""
    return (
        config.detector.width_px / 2.0
        + np.asarray(offset_mm, dtype=float)
        * 1000.0
        / config.detector.pixel_size_um
    )


def create_spectrograph_system(
    wavelength_nm: float, config: SpectrographConfig
):
    """Build the POPPY slit, lens, and detector-plane optical system."""
    o = config.optics
    system = fresnel.FresnelOpticalSystem(
        name=f"Spectrograph at {wavelength_nm} nm",
        pupil_diameter=o.pupil_diameter_mm * u.mm,
        npix=o.poppy_npix,
    )
    system.add_optic(poppy.RectangleAperture(
        width=o.slit_width_um * u.micron,
        height=o.slit_length_mm * u.mm,
        name="Long Slit",
    ))
    system.add_optic(
        poppy.QuadraticLens(
            f_lens=o.collimator_focal_length_mm * u.mm,
            name="Collimator Lens",
        ),
        distance=o.collimator_distance_mm * u.mm,
    )
    system.add_optic(
        poppy.QuadraticLens(
            f_lens=o.focusing_focal_length_mm * u.mm,
            name="Focusing Lens",
        ),
        distance=o.lens_separation_mm * u.mm,
    )
    system.add_detector(
        pixelscale=config.detector.pixel_size_um * u.micron / u.pixel,
        fov_pixels=config.detector.psf_fov_pixels,
        distance=o.detector_distance_mm * u.mm,
    )
    return system
