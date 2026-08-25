# Spectrometer-simulation

## Overview
This repository contains a numerical simulation pipeline for optical spectrographs.

## Notebooks

* **`Simulation_spec.ipynb`** - Primary simulation notebook. Simulate the output of the spectrometer configuration, as well as measuring the theoretical and simulated spectral resolution for comparison.

* **`Net_resolution.ipynb`** - Parameter sensitivity analysis. Investigate the impact of each spectrometer's parameter: slit width, diffraction grating and pixel on the overall spectral resolution (which is the full-width half maximum - FWHM). We utilize a self-calibrating spectrograph class to ensure other parameters are optimized for each fixed parameters.

> **Note:** Additional exploratory notebooks and modules are currently under development.

## Installation
Install the primary dependency, [POPPY] (Physical Optics Propagation in Python)
```bash
pip install poppy
```

## Author
Tran Vu Minh Anh

## Citations
If you use this repository or its underlying frameworks, please reference the following:
* POPPY Documentation: https://poppy-optics.readthedocs.io/
* POPPY software citation: Perrin, M. D., Soummer, R., Elliott, E. M., Lajoie, C.-P., & Sivaramakrishnan, A. (2012). POPPY: Physical Optics Propagation in Python. Astrophysics Source Code Library, ascl:1602.018.
