# Spacecraft Guidance & Navigation

> **Course assignments** — Politecnico di Milano, A.Y. 2024/25
> > Instructors: Prof. Francesco Topputto & Prof. Pierluigi Di Lizia
> > > Teaching Assistants: Sergio Bonaccorsi & Alessandro Morselli
> > >
> > > ---
> > >
> > > ## Overview
> > >
> > > This repository contains the MATLAB implementations and written reports for the two main assignments of the **Spacecraft Guidance & Navigation** course at Politecnico di Milano. The work covers two complementary macro-areas of modern astrodynamics:
> > >
> > > 1. **Guidance** — trajectory design and optimization for impulsive and low-thrust transfers, including periodic orbit computation in the Circular Restricted Three-Body Problem (CR3BP).
> > > 2. 2. **Navigation** — statistical orbit determination via batch and sequential filters, uncertainty propagation, and ground-station tracking geometry.
> > >   
> > >    3. ---
> > >   
> > >    4. ## Repository Structure
> > >   
> > >    5. ```
> > >       Spacecraft-Guidance-and-Navigation/
> > >       ├── Guidance/
> > >       │   ├── PeriodicOrbit.m          # Lagrange points & periodic orbits in CR3BP
> > >       │   ├── ImpulsiveGuidance.m      # Multi-impulse trajectory optimization (multiple shooting + Pareto)
> > >       │   ├── ContinuousGuidance.m     # Low-thrust trajectory optimization (debris-density-aware)
> > >       │   ├── kernels/                 # SPICE kernels (ex02.tm)
> > >       │   ├── mice/                    # MATLAB SPICE toolkit (MICE)
> > >       │   ├── Report.pdf               # Written report for Assignment 1
> > >       │   └── grafico_pareto_effic.png # Pareto front: ΔV vs. efficiency
> > >       │
> > >       └── Navigation/
> > >           ├── UncertaintyPropagation.m # LinCov & Unscented Transform (UT) in CR3BP
> > >           ├── BatchFilters.m           # Weighted least-squares batch orbit determination
> > >           ├── KalmanFilters.m          # EKF & UKF sequential orbit determination
> > >           ├── kernels/                 # SPICE kernels (assignment02.tm)
> > >           ├── mice/                    # MATLAB SPICE toolkit (MICE)
> > >           ├── sgp4/                    # SGP4 propagator
> > >           ├── tle/                     # Two-Line Element sets
> > >           ├── assignment02.tm          # SPICE meta-kernel
> > >           └── Report.pdf               # Written report for Assignment 2
> > >       ```
> > >
> > > ---
> > >
> > > ## Assignment 1 — Guidance
> > >
> > > ### `PeriodicOrbit.m`
> > >
> > > Computes the five **Lagrange equilibrium points** of the Earth–Moon CR3BP (mass parameter μ = 0.012150) by solving the gradient of the effective potential with `fsolve` (tolerance 1e-10). Starting from these points, the script constructs **periodic orbits** (Lyapunov / Halo families) using a differential-correction / shooting scheme that exploits the mirror symmetry of the CR3BP.
> > >
> > > ### `ImpulsiveGuidance.m`
> > >
> > > Solves an **impulsive orbital transfer** problem (≈ 1700 lines) structured around the following pipeline:
> > >
> > > - Loads SPICE ephemerides via `cspice_furnsh('ex02.tm')` and converts osculating orbital elements to Cartesian states with `par2car`.
> > > - - Applies **multiple shooting**: the transfer arc is split into sub-arcs, each propagated with `ode45`/`ode113`; interior matching conditions are enforced via `fsolve`.
> > >   - - Performs **multi-objective optimization**: a Pareto front (ΔV vs. a second performance index) is generated to expose the trade-off between propellant cost and transfer efficiency (see `grafico_pareto_effic.png`).
> > >    
> > >     - > **Note:** some solver options (`MaxIterations`, `MaxFunctionEvaluations`) are commented out to reduce runtime in the submitted version; uncommenting them reproduces the high-accuracy figures in the report's Appendix F.
> > >       >
> > >       > ### `ContinuousGuidance.m`
> > >       >
> > >       > Addresses a **low-thrust transfer** in a debris-populated environment:
> > >       >
> > >       > - Models space debris density as a radial function `q(ρ) = k₁ / (k₂ + ((ρ − ρ₀)/DU)²)` and plots it against altitude.
> > >       > - - Computes the initial circular-orbit state (position, velocity, initial mass) and propagates the **rocket equation** alongside the equations of motion.
> > >       >   - - Solves the optimal low-thrust problem via indirect methods (Euler–Lagrange / co-state equations), producing minimum-fuel or minimum-time trajectories while accounting for the debris-density penalty.
> > >       >    
> > >       >     - ---
> > >       >
> > >       > ## Assignment 2 — Navigation
> > >       >
> > >       > ### `UncertaintyPropagation.m`
> > >       >
> > >       > Compares two methods for propagating orbital uncertainty in the CR3BP:
> > >       >
> > >       > | Method | Description |
> > >       > |---|---|
> > >       > | **LinCov** | State Transition Matrix (STM) propagation: `P(t) = Φ · P₀ · Φᵀ`. Fast, accurate for small uncertainties. |
> > >       > | **UT (Unscented Transform)** | Sigma-point propagation through the full nonlinear dynamics. More accurate for larger, nonlinear uncertainties. |
> > >       >
> > >       > Initial covariance `P₀` is a 4×4 matrix (position & velocity, adimensional CR3BP units, order 1e-14 – 1e-15). A 5-point time grid between `ti` and `tf` is used to track covariance evolution.
> > >       >
> > >       > ### `BatchFilters.m`
> > >       >
> > >       > Implements **weighted least-squares batch orbit determination** using three ground stations:
> > >       >
> > >       > | Station | Location | Min. Elevation | Frequency |
> > >       > |---|---|---|---|
> > >       > | Kourou | ~0° lat | 6° | 60 passes/day |
> > >       > | Troll | Antarctica | 0° | 30 passes/day |
> > >       > | Svalbard | ~78° N | 8° | 60 passes/day |
> > >       >
> > >       > Measurement noise: σ_Az = σ_El = 125 mrad, σ_range = 10 m. SPICE is used to compute visibility windows. The filter solves `δx = (HᵀWH)⁻¹ Hᵀ W y` iterating until convergence.
> > >       >
> > >       > ### `KalmanFilters.m`
> > >       >
> > >       > Implements **sequential orbit determination** with range measurements from a lunar lander ("MOONLANDER", 78°N 15°E), tracking a satellite over a 4-hour arc (2024-11-18 16:30–20:30 UTC), sampled every 30 s (σ_ρ = 100 m):
> > >       >
> > >       > | Filter | Linearization | Dynamics |
> > >       > |---|---|---|
> > >       > | **EKF** (Extended Kalman Filter) | First-order Taylor expansion of `h(x)` | STM-based prediction |
> > >       > | **UKF** (Unscented Kalman Filter) | Sigma-point propagation (same UT as above) | Full nonlinear propagation |
> > >       >
> > >       > ---
> > >       >
> > >       > ## Dependencies & Setup
> > >       >
> > >       > | Library | Purpose | Location |
> > >       > |---|---|---|
> > >       > | [MICE (MATLAB SPICE)](https://naif.jpl.nasa.gov/naif/toolkit_MATLAB.html) | Ephemeris & frame transformations | `Guidance/mice/`, `Navigation/mice/` |
> > >       > | [SGP4](https://celestrak.org/software/tskelso-sw.php) | TLE-based orbit propagation | `Navigation/sgp4/` |
> > >       > | MATLAB R2023+ | Main runtime (requires Optimization Toolbox) | — |
> > >       >
> > >       > To run any script, ensure the `mice/src/mice/` and `mice/lib/` paths are added (the `addpath` calls at the top of each script can be uncommented if needed), and that the appropriate SPICE meta-kernel is loaded.
> > >       >
> > >       > ---
> > >       >
> > >       > ## Author
> > >       >
> > >       > **Marcello Pareschi** — MSc Aerospace Engineering, Politecnico di Milano
> > >       > [GitHub @PaneeVino](https://github.com/PaneeVino)
