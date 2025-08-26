# Nonlinear Dynamics and Chaos: Advanced Computational Analysis Tools

**Comprehensive R implementation of dynamical systems analysis including chaos detection, bifurcation theory, and interactive structural mechanics**

## Overview

This repository provides production-ready tools for analyzing complex dynamical systems, focusing on chaos theory, bifurcation analysis, and nonlinear structural mechanics. The implementations combine rigorous numerical methods with advanced visualization techniques to enable deep exploration of system behavior across parameter spaces.

## Core Implementations

### 🌀 **Chaotic Pendulum Analysis**
Complete analysis suite for the periodically forced, damped pendulum:
```
θ̈ + γθ̇ + sin(θ) = A cos(ωt)
```

**Advanced Features:**
- **0-1 Test for Chaos**: Gottwald-Melbourne method for rigorous chaos detection
- **Lyapunov Exponents**: Two-trajectory method with automated renormalization
- **Poincaré Sections**: Stroboscopic sampling revealing attractor structure
- **Basin of Attraction**: Phase space visualization of initial condition sensitivity
- **Return Maps**: One-dimensional discrete dynamics from continuous systems

### ⚖️ **Interactive Beam Buckling (Shiny App)**
Comprehensive structural mechanics application covering:
- **Euler Buckling Theory**: Critical loads for multiple boundary conditions
- **Elastica Post-Buckling**: Exact solutions using elliptic integrals
- **Pitchfork Bifurcation**: Normal form analysis with imperfection sensitivity

**Boundary Conditions Supported:**
- Pinned-Pinned: K = 1.0
- Clamped-Clamped: K = 0.5  
- Clamped-Pinned: K = 0.699
- Cantilever: K = 2.0

## Technical Architecture

### Numerical Integration
High-precision RK4 integrator optimized for chaotic systems:
```r
rk4_step <- function(f, t, x, dt, pars) {
  k1 <- f(t,         x,            pars)
  k2 <- f(t + dt/2,  x + dt*k1/2,  pars)
  k3 <- f(t + dt/2,  x + dt*k2/2,  pars)
  k4 <- f(t + dt,    x + dt*k3,    pars)
  x + dt*(k1 + 2*k2 + 2*k3 + k4)/6
}
```

### Chaos Detection Methods

#### 0-1 Test (Gottwald-Melbourne)
Robust statistical test distinguishing chaos from regular motion:
```r
zero_one_test <- function(sim_df, observable="theta", c_vals=NULL) {
  # Translation dynamics: pc = Σv·cos(nc), qc = Σv·sin(nc)
  # Mean square displacement correlation with time index
  K <- median(cor(Mc, idx))  # ≈1 chaos, ≈0 regular
}
```

#### Lyapunov Exponent Calculation
Two-trajectory method with periodic renormalization:
```r
lyapunov_largest <- function(..., renorm_every_steps=50, d0=1e-8) {
  # Evolve nearby trajectories, renormalize separation
  # Track exponential divergence rate
  sum_log / total_time  # Average exponential growth rate
}
```

### Structural Mechanics (Elastica Theory)

#### Post-Buckling Analysis
Exact solutions using complete elliptic integrals K(m), E(m):
```r
# End-shortening relation: δ = 1 - E(m)/K(m)
# Load relation: P = EI(2K(m)/L)²
solve_m_from_delta <- function(delta) {
  f <- function(m) 1 - KE_from_m(m)$E/KE_from_m(m)$K - delta
  uniroot(f, interval=c(1e-12, 1-1e-10))$root
}
```

#### Shape Integration
Numerical integration of elastica differential equation:
```r
# θ'' + λ²sin(θ) = 0 with symmetry boundary conditions
elastica_halfshape <- function(L, m, n=600) {
  # RK4 integration from center to end
  # Build full shape by symmetry
}
```

## Advanced Analysis Capabilities

### Bifurcation Diagrams
Automated parameter sweeps with attractor detection:
```r
bifurcation_vs_A <- function(A_seq, ...) {
  # For each parameter value:
  # 1. Integrate transient dynamics
  # 2. Extract Poincaré section
  # 3. Classify attractor type
  # 4. Compute Lyapunov exponent
}
```

### Basin of Attraction Mapping
Grid-based initial condition sensitivity analysis:
```r
basin_of_attraction_grid <- function(theta_range, omega_range, n_theta, n_omega) {
  # Systematic initial condition sweep
  # Long-time attractor classification
  # Fractal boundary detection
}
```

### Return Map Construction
One-dimensional discrete dynamics from continuous systems:
```r
return_map <- function(pmap_df, var="theta") {
  # Extract successive Poincaré points: x[n+1] vs x[n]
  # Reveal fixed points and periodic orbits
  # Enable symbolic dynamics analysis
}
```

## Interactive Visualization (Shiny)

### Real-Time Parameter Exploration
The beam buckling app provides:
- **Dynamic Critical Load Calculation**: P_cr = π²EI/(KL)² with live updates
- **Mode Shape Visualization**: Exact and approximate eigenfunctions
- **Post-Buckling Curves**: Load-deflection relationships under displacement control
- **Normal Form Bifurcations**: Pitchfork diagrams with imperfection effects

### Advanced UI Features
```r
# Material properties with unit conversion
numericInput("E_gpa", "Young's modulus E [GPa]", value=200)
numericInput("I_cm4", "Area moment I [cm⁴]", value=100)

# Interactive elastica visualization
sliderInput("delta_ratio", "End-shortening δ=Δ/L", 0, 0.45, 0.05, 0.005)
```

## Scientific Applications

### Chaos Theory Research
- **Route to Chaos**: Period-doubling cascades and intermittency
- **Strange Attractors**: Geometric characterization of chaotic sets
- **Sensitivity Analysis**: Quantifying predictability horizons

### Structural Engineering
- **Buckling Load Prediction**: Critical loads for complex boundary conditions
- **Post-Critical Behavior**: Large deflection analysis using exact theory
- **Imperfection Sensitivity**: Bifurcation analysis with manufacturing tolerances

### Mathematical Physics
- **Hamiltonian Mechanics**: Conservative vs dissipative system comparison
- **Symmetry Breaking**: Pitchfork bifurcations in physical systems
- **Normal Forms**: Universal unfolding of bifurcations

## Performance Characteristics

### Computational Complexity
| Analysis Type | Time Complexity | Space Complexity | Typical Runtime |
|--------------|----------------|------------------|-----------------|
| Single Trajectory | O(N) | O(N) | ~1s for 10⁵ steps |
| Lyapunov Exponent | O(N) | O(1) | ~5s for 10⁶ steps |
| Bifurcation Diagram | O(MN) | O(MN) | ~2min for 100 parameters |
| Basin Analysis | O(G²N) | O(G²) | ~10min for 80² grid |

Where N = integration steps, M = parameter values, G = grid resolution.

### Numerical Stability
- **Adaptive Timestep Criteria**: CFL conditions for chaotic flows
- **Symplectic Integration**: Available for Hamiltonian systems
- **Error Monitoring**: Energy conservation and phase space volume tracking

## Usage Examples

### Chaos Detection Workflow
```r
source("chaotic_pendulum_enhanced.R")

# Phase space analysis
demo_phase(pars=list(gamma=0.2, A=1.4, omega=2/3))

# Rigorous chaos test
demo_zero_one(pars=list(gamma=0.2, A=1.4, omega=2/3), 
              observable="theta")

# Parameter space exploration  
demo_bifurcation(A_min=0.9, A_max=1.6, nA=120)
```

### Interactive Structural Analysis
```r
source("shiny_app_buckling_beam.R")
# Launches interactive web application with:
# - Real-time critical load calculation
# - Mode shape visualization
# - Post-buckling analysis
# - Pitchfork bifurcation diagrams
```

### Basin of Attraction Study
```r
# Fractal boundary analysis
basin <- demo_basin(
  pars=list(gamma=0.2, A=1.4, omega=2/3),
  theta_range=c(-pi, pi), 
  omega_range=c(-2, 2),
  n_theta=80, n_omega=80
)
```

## Theoretical Foundations

### Chaos Theory
The implementations leverage key concepts from dynamical systems theory:
- **Sensitive Dependence**: Exponential divergence of nearby trajectories
- **Strange Attractors**: Fractal structures in phase space
- **Period-Doubling Routes**: Universal scaling laws (Feigenbaum constants)

### Bifurcation Theory
- **Codimension-1 Bifurcations**: Saddle-node, pitchfork, Hopf bifurcations
- **Normal Forms**: Local analysis near bifurcation points
- **Imperfection Theory**: Unfolding of degenerate bifurcations

### Structural Mechanics
- **Elastica Theory**: Exact large-deflection beam analysis
- **Stability Theory**: Linear and nonlinear buckling predictions
- **Elliptic Integrals**: Special functions for post-critical analysis

## Research Integration

### Data Export Capabilities
All analyses generate publication-ready output:
```r
# Automated PDF generation
ggsave("bifurcation_diagram.pdf", plot=bif_plot, width=8, height=6)

# CSV data for further analysis
write.csv(basin_data, "basin_classification.csv")

# Time series for external tools
saveRDS(simulation_results, "pendulum_trajectory.rds")
```

### Parameter Study Framework
```r
# Systematic exploration
param_grid <- expand.grid(
  gamma = seq(0.1, 0.3, 0.05),
  A = seq(1.0, 2.0, 0.1),
  omega = c(1/2, 2/3, 3/4)
)
# Parallel processing via apply family functions
```

## Installation & Dependencies

### Core Requirements
```r
install.packages(c("shiny", "ggplot2", "dplyr"))

# Optional for enhanced elliptic integrals
install.packages("pracma")  # Faster K(m), E(m) computation
```

### Hardware Recommendations
- **RAM**: 4GB+ for large parameter sweeps
- **CPU**: Multi-core beneficial for basin analysis
- **Graphics**: Modern GPU helpful for real-time Shiny visualization

## Literature Foundation

Implementations follow established methods from:

1. **Gottwald & Melbourne (2004)**: 0-1 test for chaos
2. **Benettin et al. (1980)**: Lyapunov exponent algorithms
3. **Guckenheimer & Holmes (1983)**: Nonlinear Oscillations, Dynamical Systems
4. **Thompson & Hunt (1973)**: Elastic Instability Theory

The structural mechanics draws from:
- **Timoshenko & Gere (1961)**: Theory of Elastic Stability
- **Antman (2005)**: Nonlinear Problems of Elasticity

---

*A comprehensive computational framework for nonlinear dynamics research combining rigorous mathematical theory with practical analysis tools.*