<div align="center">
  <img src="pics/logo_light.svg#gh-light-mode-only" alt="WannierBSE Logo" width="700" height="auto">
  <img src="pics/logo_dark.svg#gh-dark-mode-only" alt="WannierBSE Logo" width="700" height="auto">

# WannierBSE: High-Performance Bethe-Salpeter Equation Solver

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-Supported-blue.svg)](https://www.mathworks.com/products/MATlab.html)
[![GitHub Stars](https://img.shields.io/github/stars/WannierBSE/WannierBSE-MATLAB-?style=social)](https://github.com/WannierBSE/WannierBSE-MATLAB-)

[**What is WannierBSE?**](#what-is-wannierbse) | [**Getting Started**](#getting-started) | [**Technical Workflow**](#technical-workflow) | [**Examples**](#examples) | [**Citation**](#citation) | [**Credits**](#credits) | [**Supported By**](#supported-by)
</div>

---

## What is WannierBSE?

**WannierBSE** is a high-performance computational suite engineered for solving the **Bethe-Salpeter Equation (BSE)** within the **Wannier Tight-Binding (WTB)** framework. Specifically optimized for 2D materials, WannierBSE bridges the gap between first-principles accuracy and model-scale efficiency.

By leveraging MATLAB’s robust numerical engine, the package enables researchers to simulate excitonic spectra, dipoles, and light-matter interactions under complex dielectric environments.

## Getting Started

To get started with WannierBSE, please refer to our comprehensive documentation:

*   [**User Guide (PDF)**](WannierBSE_v1.0_User_Guide.pdf) - Detailed technical manual and theory.
*   [**QuickStart Guide (PDF)**](WannierBSE_v1.0_QuickStart_Guide.pdf) - Brief setup and first run instructions.
*   [**Official Website**](https://quantum.web.nycu.edu.tw/wannierbse) - News and updates.

### Prerequisites
*   **MATLAB** (Recommended: R2020b or later)
*   **Parallel Computing Toolbox** (Highly recommended for large-scale calculations)

## Technical Workflow

WannierBSE features a modular data pipeline designed to handle both raw generation and pre-computed data ingestion.

### 1. Data Preparation
Users can initiate calculations through two primary pathways:

> 🔄 **Generation Route:** Place `wannier90_hr.dat` in the `/User_input/` folder. The solver will automatically generate tight-binding bands, k-meshes, and dielectric functions, saving them in the `/Precomputed_data/` folder.
> 
> 📥 **Custom/Direct Route:** Import user-defined files (e.g., `v*_TB.txt`, `c*_TB.txt`, `kmesh.txt`) directly into `/User_input/` to bypass internal generation.

### 2. Configuration & Control
The simulation environment is governed by dedicated control files located in the `/Parameters/` directory:
*   `control.txt` & `structure.txt`
*   `WF_centers.txt` & `kmesh_control.txt`
*   `WTB_control.txt` & `dielectric_control.txt`

### 3. Execution & Visualization
1.  **Core Solver:** Execute `WBSE.m` in MATLAB. The script utilizes an optimized, symmetry-aware Hamiltonian constructor and a high-performance parallelized solver.
2.  **Post-Processing:** Use the auxiliary script `EX_plot.m` to visualize the exciton energy spectrum.
3.  **Data Management:** All results, including energy spectra (`Ex.mat`) and wavefunctions (`A.mat`), are automatically saved in the `/Exciton_data/` folder.

## Examples

To help you get started quickly, we provide several pre-configured examples in the repository. These include sample input files and parameters for typical 2D materials.

📂 [**Browse Examples in the Repository**](https://github.com/WannierBSE/WannierBSE-MATLAB-/tree/main/Examples)

Each example includes:
*   User input and precomputed data files
*   Ready-to-use control parameter files
*   Expected results.

## Citation

To support the continued development of WannierBSE, please cite the following references:

> **Primary Reference:**
> Peng, G.-H. et al. "Distinctive Signatures of the Spin- and Momentum-Forbidden Dark Exciton States in the Photoluminescence of Strained WSe2 Monolayers under Thermalization", 
> *Nano Lett.* 19, 4, 2299–2312 (2019). 
> [https://doi.org/10.1021/acs.nanolett.8b04786](https://doi.org/10.1021/acs.nanolett.8b04786)

> **Non-local Dielectric Function:**
> Li, W.-H. et al. "The Key Role of Non-Local Screening in the Environment-Insensitive Exciton Fine Structures of Transition-Metal Dichalcogenide Monolayers", 
> *Nanomaterials* 13, 11, 1739 (2023). 
> [https://doi.org/10.3390/nano13111739](https://doi.org/10.3390/nano13111739)

## Credits

Developed by the research group of **Prof. Shun-Jen Cheng** at National Yang Ming Chiao Tung University (NYCU).

**Core Developers:**
*   Dr. Ping-Yuan Lo, Dr. Wei-Hua Li, Dr. Guan-Hao Peng, Dr. Jhen-Dong Lin, Dr. Vo Khuong Dien, Dr. Oscar Javier Gomez Sanchez, Mr. Ching-Hung Shih, Mr. Kun-Yi Lin, and Prof. Shun-Jen Cheng.

## Supported By

We gratefully acknowledge the support from the following organizations:

<div align="center">
  <img src="pics/NSTC-Light.webp#gh-light-mode-only" height="40" align="left">
  <img src="pics/NSTC-Dark.webp#gh-dark-mode-only" height="40" align="left">
  
  <img src="pics/TRIAD%20LIGHT%20INNOVATION-Light.webp#gh-light-mode-only" height="45" align="right">
  <img src="pics/TRIAD%20LIGHT%20INNOVATION-Dark.webp#gh-dark-mode-only" height="45" align="right">
  
  <img src="pics/NCHC-Light.webp#gh-light-mode-only" height="40">
  <img src="pics/NCHC-Dark.webp#gh-dark-mode-only" height="40">
</div>
<br clear="all">

---
<div align="center">
  Copyright © 2026 The WannierBSE Development Team. All rights reserved.
</div>
