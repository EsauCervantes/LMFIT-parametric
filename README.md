# LMFIT Parametric Scan for Higgs-Portal Dark Matter

This repository contains Python scripts used to perform numerical parameter scans for the particle-physics model studied in:

**E. Cervantes, O. Perez-Figueroa, R. Perez-Martinez, S. Ramos-Sanchez,  
"Higgs-portal dark matter from non-supersymmetric strings",  
Phys. Rev. D 107, 115007 (2023), arXiv:2302.08520.**

The code was developed as part of the phenomenological analysis of Higgs-portal dark matter candidates arising from non-supersymmetric string-inspired models. It combines spectrum generation, vacuum-stability checks, relic-density calculations, and numerical minimization/scanning of the model parameter space.

## Purpose

The scripts implement a workflow for scanning a 14-dimensional parameter space and comparing generated model points against selected phenomenological constraints, including:

- the observed Higgs mass,
- the electroweak vacuum expectation value,
- the dark matter relic abundance,
- tree-level unitarity constraints,
- vacuum stability or metastability,
- consistency of the generated vacuum with the input value of `tan(beta)`.

The objective function is expressed as a residual vector, which is minimized using routines from [`lmfit`](https://lmfit.github.io/lmfit-py/). Some versions of the scan also include exploratory Monte Carlo / MCMC-style scans using `emcee`.

## External tools

The workflow interfaces with the following external particle-physics codes:

- **SPheno** — spectrum generation,
- **Vevacious / VevaciousPlusPlus** — vacuum-stability analysis,
- **MicrOmegas** — dark matter relic-density calculation.

These tools are not included in this repository and must be installed separately.

## Repository structure

The main scripts are:

- `scan_chisqrd.0.1.py`  
  Main driver script for the chi-square minimization / parameter scan.

- `spvevmicro.py`  
  Helper module that prepares input files, runs SPheno, Vevacious, and MicrOmegas, parses their outputs, and computes the residual vector or chi-square value.

- `spvevmicro_tanbeta.py`  
  Variant of the helper module including an additional consistency contribution associated with `tan(beta)`.

- `chisqrd_tanbeta/`, `exp_step/`, `paso2/`, `old/`  
  Development folders containing alternative versions, exploratory scans, or older scan stages.

## Important note

This repository is intended mainly as a research-code archive accompanying the publication above. The scripts contain hard-coded paths and were written for the local software environment used during the original analysis. They may require modification before they can be run on a different machine.

In particular, users should update the paths to:

- the SPheno executable and input files,
- the Vevacious executable and XML input files,
- the MicrOmegas model directory and executable,
- the corresponding output files used by the parser.

## Citation

If you use or refer to this code, please cite:

```bibtex
@article{Cervantes:2023,
  author = {Cervantes, Esau and Perez-Figueroa, Omar and Perez-Martinez, Ricardo and Ramos-Sanchez, Saul},
  title = {Higgs-portal dark matter from non-supersymmetric strings},
  journal = {Phys. Rev. D},
  volume = {107},
  pages = {115007},
  year = {2023},
  eprint = {2302.08520},
  archivePrefix = {arXiv},
  primaryClass = {hep-ph}
}
