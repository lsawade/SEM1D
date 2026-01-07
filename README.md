# 1D Wave Equation - Spectral Element Method

MATLAB implementation of 1D wave propagation using the Spectral Element Method (SEM) with real-time visualization.

## Directory Structure

```
SEM1D/
├── wave.m                    # Main simulation with real-time visualization
├── README.md                 # This file
├── output/                   # Snapshot files directory
└── functions/                # Required library functions (12 files)
    ├── define_derivative_matrix.m
    ├── lagrange_deriv_GLL.m
    ├── zwgljd.m, zwgjd.m
    ├── jacg.m, jacobf.m
    ├── pnleg.m, pndleg.m
    ├── pnormj.m, gammaf.m
    └── endw1.m, endw2.m
```

> **Note:** The `functions/` directory is automatically added to MATLAB's path when `wave.m` is executed. No manual path management is required.

## Usage

### Running the Simulation

**Dirichlet Boundary Conditions (default):**
```matlab
wave()
% or explicitly
wave('dirichlet')
```

**Neumann Boundary Conditions:**
```matlab
wave('neumann')
```

The simulation displays real-time animation showing wave propagation. Snapshots are saved to `output/` every 25 timesteps.

## Solution Implementation

### 1. Global Mass Matrix
Assembles global mass matrix from local element contributions:
```matlab
mass_local = wgll(i) * rho(i, ispec) * jacobian(i, ispec);
mass_global(iglob) = mass_global(iglob) + mass_local;
```

### 2. Time Integration (Newmark Scheme)

**Time step parameters:**
```matlab
dt = 0.25;
halfdt = 0.5 * dt;
halfdt2 = 0.5 * dt * dt;
```

**Predictor step:**
```matlab
displ = displ + dt * veloc + halfdt2 * accel;
veloc = veloc + halfdt * accel;
accel = zeros(NGLOB, 1);
```

**Corrector step:**
```matlab
accel = accel + daccel;
veloc = veloc + halfdt * accel;
```

### 3. Initial Conditions
Gaussian pulse centered at x = 50:
```matlab
displ = exp(-0.1 * (x - 50.0).^2);
veloc = zeros(NGLOB, 1);
accel = zeros(NGLOB, 1);
```

### 4. Boundary Conditions
- **Dirichlet:** Fixes displacement at both ends to zero
- **Neumann:** Traction-free (natural) boundary conditions

### 5. Solver
Assembles global force vector and solves for acceleration:
```matlab
stiffness_local = -sum(wgll .* shear(:, ispec) .* ...
                     inverse_jacobian(:, ispec) .* ...
                     hprime(i, :)' .* hprime(j, :)');
daccel = F_global ./ mass_global;
```

## Heterogeneous Models

The code supports **spatially variable material properties** through user-defined functions. By default, the model is homogeneous with ρ = 1 and μ = 1.

### Defining Material Properties

Edit the material property functions in `wave.m` (lines 70-73):

```matlab
% Density as a function of x
get_density = @(x) ones(size(x));  % Default: homogeneous

% Shear modulus as a function of x
get_shear = @(x) ones(size(x));    % Default: homogeneous
```

### Example Models

**Two-layer model** (interface at x = 50):
```matlab
get_density = @(x) 1.0 + 0.5 * (x > 50);  % ρ₁=1.0, ρ₂=1.5
get_shear = @(x) 1.0 + 1.0 * (x > 50);    % μ₁=1.0, μ₂=2.0
```

**Gradient model:**
```matlab
get_density = @(x) 1.0 + 0.01 * x;
get_shear = @(x) 1.0 + 0.02 * x;
```

**Low-velocity zone:**
```matlab
get_density = @(x) 1.0 - 0.3 * exp(-0.01 * (x - 50).^2);
get_shear = @(x) 1.0 - 0.5 * exp(-0.01 * (x - 50).^2);
```

A figure showing the material properties will be displayed when the simulation starts.

## Parameters

Simulation parameters defined in `wave.m`:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `NSPEC` | 20 | Number of spectral elements |
| `NGLL` | 7 | Number of GLL points per element |
| `NGLOB` | 121 | Total number of global points |
| `NSTEP` | 2000 | Number of time steps |
| `DT` | 0.25 | Time step (seconds) |
| `LENGTH` | 100 | Domain length |

**Total simulation time:** 500 seconds (2000 steps × 0.25 s)

**Material properties** (ρ, μ) are defined via functions `get_density(x)` and `get_shear(x)`. See [Heterogeneous Models](#heterogeneous-models) section.
