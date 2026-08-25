# Spectrometer-simulation

## Overview
This repository contains a numerical simulation pipeline for optical spectrographs.

## Notebooks

* **`Simulation_spec.ipynb`** - This notebook is the primary notebook to simulate the output of the spectrometer configuration, as well as measuring the theoretical and simulated spectral resolution for comparison.

* **`Net_resolution.ipynb`** - This notebook is used to investigate the impact of each spectrometer's parameter: slit width, diffraction grating and pixel on the overall spectral resolution (which is the full-width half maximum - FWHM). We utilize a self-calibrating spectrograph class to ensure other parameters are optimized for each fixed parameters.

> **Note:** Additional exploratory notebooks and modules are currently under development.

## Installation
Installation of the POPPY library is necessary:
```
pip install poppy
```

## Author
Tran Vu Minh Anh

## Citations
https://poppy-optics.readthedocs.io/
