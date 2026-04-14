# Spacecraft Guidance & Navigation

A collection of MATLAB implementations and written reports for the **Spacecraft Guidance & Navigation** course at [Politecnico di Milano](https://www.polimi.it), A.Y. 2024/25.

> **Instructors:** Prof. Francesco Topputto & Prof. Pierluigi Di Lizia
> **Teaching Assistants:** Sergio Bonaccorsi & Alessandro Morselli

---

## Overview

This repository contains the course assignments covering two complementary macro-areas of modern astrodynamics:

| Macro-area | Topics |
|---|---|
| **Guidance** | Trajectory design and optimization for impulsive and low-thrust transfers, periodic orbit computation in the CR3BP |
| **Navigation** | Statistical orbit determination via batch and sequential filters, uncertainty propagation, ground-station tracking geometry |

Both assignments include a MATLAB implementation and a written report.

---

## Repository Structure

```
Spacecraft-Guidance-and-Navigation/
├── Guidance/
│   ├── PeriodicOrbit.m          # Lagrange points & periodic orbits in CR3BP
│   ├── ImpulsiveGuidance.m      # Multi-impulse trajectory optimization (multiple shooting + Pareto)
│   ├── ContinuousGuidance.m     # Low-thrust trajectory optimization (debris-density-aware)
│   ├── kernels/                 # SPICE kernels (ex02.tm)
│   ├── mice/                    # MATLAB SPICE toolkit (MICE)
│   ├── Report.pdf               # Written report for Assignment 1
│   └── grafico_pareto_effic.png # Pareto front: ΔV vs. efficiency
└── Navigation/
    ├── UncertaintyPropagation.m # LinCov & Unscented Transform (UT) in CR3BP
    ├── BatchFilters.m           # Weighted least-squares batch orbit determination
    ├── KalmanFilters.m          # EKF & UKF sequential orbit determination
    ├── kernels/                 # SPICE kernels (assignment02.tm)
    ├── mice/                    # MATLAB SPICE toolkit (MICE)
    ├── sgp4/                    # SGP4 propagator
    ├── tle/                     # Two-Line Element sets
    ├── assignment02.tm          # SPICE meta-kernel
    └── Report.pdf               # Written report for Assignment 2
```

---

## Features

### Assignment 1 — Guidance

#### PeriodicOrbit.m

Computes the five Lagrange equilibrium points of the Earth–Moon CR3BP (mass parameter μ = 0.012150) by solving the gradient of the effective potential with `fsolve` (tolerance 1e-10). Starting from these points, the script constructs periodic orbits (Lyapunov / Halo families) using a differential-correction / shooting scheme that exploits the mirror symmetry of the CR3BP.

#### ImpulsiveGuidance.m

Solves an impulsive orbital transfer problem structured around the following pipeline:

- **SPICE ephemerides** loaded via `cspice_furnsh('ex02.tm')`; osculating orbital elements converted to Cartesian states with `par2car`.
- **Multiple shooting**: the transfer arc is split into sub-arcs, each propagated with `ode45`/`ode113`; interior matching conditions enforced via `fsolve`.
- **Multi-objective optimization**: a Pareto front (ΔV vs. a second performance index) exposes the trade-off between propellant cost and transfer efficiency (see `grafico_pareto_effic.png`).

#### ContinuousGuidance.m

Addresses a low-thrust transfer in a debris-populated environment:

- **Debris density model**: radial function q(ρ) = k₁ / (k₂ + ((ρ − ρ₀)/DU)²).
- **Rocket equation** propagated alongside the equations of motion.
- **Indirect optimal control** (Euler–Lagrange / co-state equations) producing minimum-fuel or minimum-time trajectories with a debris-density penalty.

---

### Assignment 2 — Navigation

#### UncertaintyPropagation.m

Compares two methods for propagating orbital uncertainty in the CR3BP:

| Method | Description |
|---|---|
| **LinCov** | State Transition Matrix (STM) propagation: P(t) = Φ · P₀ · Φᵀ. Fast, accurate for small uncertainties. |
| **UT (Unscented Transform)** | Sigma-point propagation through the full nonlinear dynamics. More accurate for larger, nonlinear uncertainties. |

Initial covariance P₀ is a 4×4 matrix (position & velocity, adimensional CR3BP units, order 1e-14 – 1e-15). A 5-point time grid between tᵢ and tᶠ is used to track covariance evolution.

#### BatchFilters.m

Implements weighted least-squares batch orbit determination using three ground stations:

| Station | Location | Min. Elevation | Frequency |
|---|---|---|---|
| Kourou | ~0° lat | 6° | 60 passes/day |
| Troll | Antarctica | 0° | 30 passes/day |
| Svalbard | ~78° N | 8° | 60 passes/day |

Measurement noise: σ_Az = σ_El = 125 mrad, σ_range = 10 m. SPICE is used to compute visibility windows. The filter solves δx = (HᵀWH)⁻¹ Hᵀ W y iterating until convergence.

#### KalmanFilters.m

Implements sequential orbit determination with range measurements from a lunar lander ("MOONLANDER", 78°N 15°E), tracking a satellite over a 4-hour arc (2024-11-18 16:30–20:30 UTC), sampled every 30 s (σ_ρ = 100 m):

| Filter | Linearization | Dynamics |
|---|---|---|
| **EKF** (Extended Kalman Filter) | First-order Taylor expansion of h(x) | STM-based prediction |
| **UKF** (Unscented Kalman Filter) | Sigma-point propagation (same UT as above) | Full nonlinear propagation |

---

## Getting Started

### Requirements

- MATLAB R2023b or later (requires Optimization Toolbox).
- MICE (MATLAB SPICE Toolkit) — included in `Guidance/mice/` and `Navigation/mice/`.
- SGP4 propagator — included in `Navigation/sgp4/`.

### Running an Example

1. Clone the repository:
```
git clone https://github.com/PaneeVino/Spacecraft-Guidance-and-Navigation.git
```

2. Open MATLAB and add the MICE paths (each script contains `addpath` calls at the top that can be uncommented):
```matlab
addpath(genpath('Guidance/mice/src/mice/'))
addpath(genpath('Guidance/mice/lib/'))
```

3. Load the appropriate SPICE meta-kernel and run the desired script:
```matlab
cspice_furnsh('Guidance/kernels/ex02.tm')
run('Guidance/PeriodicOrbit.m')
```

---

## Academic Context

These codes were developed as part of the coursework for the **Spacecraft Guidance & Navigation** course at Politecnico di Milano. Benchmark problems include classical astrodynamics test cases covering Halo orbit continuation in the CR3BP, multi-impulse transfer optimisation, and statistical orbit determination with realistic ground-station geometries.

---

## License

This repository is shared for educational purposes. If you use or adapt this code, please acknowledge the source.
