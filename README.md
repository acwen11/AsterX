<img align="top" src="Docs/figures/asterx.png" width="140">

**AsterX** is a GPU-accelerated GRMHD code for dynamical spacetimes, written in C++. It is built on the [CarpetX](https://github.com/eschnett/CarpetX) driver, which is designed for use with the [Einstein Toolkit](https://einsteintoolkit.org/). **CarpetX** itself is based on [AMReX](https://amrex-codes.github.io), a software framework for block-structured adaptive mesh refinement (AMR).

[![GitHub CI](https://github.com/jaykalinani/AsterX/workflows/CI/badge.svg)](https://github.com/jaykalinani/AsterX/actions)  
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://github.com/jaykalinani/AsterX/blob/main/LICENSE.md)

## Overview

* Heavily derived from the GRMHD code [Spritz](https://zenodo.org/records/15440856)  
* Solves the GRMHD equations in 3D Cartesian coordinates on dynamical spacetimes using high-resolution shock-capturing (HRSC) schemes  
* Based on the flux-conservative Valencia formulation  
* Directly evolves the staggered vector potential  

## Available Modules

* `AsterX` – Core GRMHD module  
* `AsterMasks` – Provides masking functionality  
* `AsterSeeds` – Initial data module  
* `AsterUtils` – Utility functions  
* `Con2PrimFactory` – Conservative-to-primitive variable recovery routines  
* `EOSX` – Equation of state driver  
* `FishboneMoncriefIDX` – Initial data for Fishbone–Moncrief disks  
* `ID_TabEOS_HydroQuantities` – Initializes hydrodynamic quantities for tabulated EOS  
* `ReconX` – Reconstruction scheme module  
* `TOVSolverX` – Modified version of the TOVSolver thorn from the Einstein Toolkit  

## Getting Started

* Instructions for downloading and building AsterX with the Einstein Toolkit are available [here](https://github.com/EinsteinToolkit/CarpetX/wiki/Getting-Started).  
* Simfactory files for various clusters and setup instructions can be found [here](https://github.com/lwJi/ETK-Compile-Guides).  
* Example Jupyter notebooks and plotting scripts are available [here](https://github.com/jaykalinani/AsterX-Docs/tree/ETX_2024_06).  

## Useful Repositories

* [CarpetX](https://github.com/EinsteinToolkit/CarpetX) – Next-generation driver for the Einstein Toolkit  
* [SpacetimeX](https://github.com/EinsteinToolkit/SpacetimeX) – Modules for spacetime evolution 
* [BNSTools](https://github.com/jaykalinani/BNSTools) – Utilities supporting BNS merger simulations  
* [nuX](https://github.com/jaykalinani/nuX) – Modules for an upcoming neutrino transport code  

## Code Papers

* [Kalinani et al. (2024)](https://iopscience.iop.org/article/10.1088/1361-6382/ad9c11) – Introduction to the AsterX code  
* [Ji et al. (2025)](https://arxiv.org/abs/2503.09629) – Subcycling algorithm and spacetime solver in CarpetX  

## Related Talks and Tutorials

* [Using CarpetX: A Guide for Early Adopters](http://einsteintoolkit.org/seminars/2021_03_18/index.html)  
  Seminar by Erik Schnetter providing an overview of **CarpetX** capabilities.

* [Tutorial: GPUs and the Einstein Toolkit](https://einsteintoolkit.github.io/et2022uidaho/lectures/38-Tutorial8/index.html)  
  Tutorial by Lorenzo Ennoggi, Jay Kalinani, and Federico Lopez Armengol during the North American Einstein Toolkit Workshop 2022. Includes an introduction to **AsterX** and hands-on sessions.

* [AsterX: A New Open-Source GPU-Accelerated GRMHD Code for Dynamical Spacetimes](https://drive.google.com/file/d/1Z4i--W56mxeNIu598LQTpEEowX56FOoD/view?usp=sharing)  
  Slides from a talk by Jay Kalinani at the APS April Meeting 2023.

