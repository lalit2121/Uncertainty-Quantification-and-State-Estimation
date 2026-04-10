##
![Screenshot](Screenshot(54).png)

# Uncertainty Quantification and State Estimation in the Evolution of Reachable Sets of Maneuvers for Space Mission Applications

**Università di Pisa**  
**Dipartimento di Ingegneria Civile e Industriale**  
**Laurea Magistrale in Ingegneria Aerospaziale**

**Candidate:** Lalit Deshmukh  
**Advisors:** Prof. Gianpietro Di Rito (University of Pisa) & Prof. Simone Servadio (Iowa State University)

---

## Overview

This repository contains the complete implementation and supporting materials for a novel computational framework that integrates **Differential Algebra (DA)** with **Multiple Model Adaptive Estimation (MMAE)** to perform robust orbit determination and maneuver detection in the presence of unknown impulsive velocity changes.

The framework addresses a critical challenge in modern space operations: accurately estimating spacecraft states and reachable sets when non-cooperative objects execute sudden, unannounced maneuvers. Traditional single-model filters fail under these conditions. By combining DA’s high-order polynomial uncertainty propagation with a Bayesian multi-model adaptive filter, the system achieves fast, deterministic uncertainty quantification and real-time maneuver identification—capabilities essential for space situational awareness, collision avoidance, and autonomous navigation.

The code is implemented in C++ and leverages the DACE library for differential-algebraic computations, an RK78 integrator for high-fidelity orbital propagation, and a bank of parallel extended Kalman filters for adaptive estimation.

---

## Key Features

- **DA-based uncertainty propagation** – Second-order Taylor expansions propagate nonlinear orbital dynamics and uncertainties analytically, achieving ~3× faster computation than Monte Carlo methods while maintaining comparable accuracy for moderate perturbations.
- **MMAE for unknown maneuvers** – A bank of 614 parallel models hypothesizes possible velocity impulses (0.1–0.3 km/s in a discretized spherical grid). Bayesian weight updates and pruning dynamically focus on the most probable maneuver scenarios.
- **Real-time orbit determination filter** – DA-enhanced EKF handles nonlinear measurements (range, azimuth, elevation) without linearization approximations.
- **Monte Carlo validation suite** – 100–1000 independent runs quantify statistical performance, error bounds, and 3σ confidence regions.
- **Reachable-set visualization** – 2D (XY) and 3D projections of post-maneuver uncertainty envelopes.
- **Model success metrics** – Automatic classification of detection performance (perfect match, partial match, opposite-direction, or failure) with overall success rate of ~91.1 %.

---

## Repository Structure
## Repository Structure
├── src/
│   ├── DA_Orbital_Perturbations/      # Chapter 3: DA vs. Monte Carlo comparison
│   ├── DA_EKF_Orbit_Determination/    # Baseline DA-enhanced EKF filter
│   ├── MMAE_Impulse_Estimation/       # Core multi-model adaptive filter (614 models)
│   └── utils/                         # RK78 integrator, measurement models, CSV output
├── results/
│   ├── figures/                       # All plots from the thesis (3D reachable sets, error vs. time, model weights, etc.)
│   ├── csv/                           # Raw Monte Carlo statistics, model probabilities, 3σ bounds
│   └── success_analysis/              # Automatic success-rate classification
├── docs/
│   ├── thesis_abstract.md
│   └── methodology_notes.md
├── CMakeLists.txt
├── README.md
└── LICENSE


