MATLAB Version of HW8A - 1D Wave Equation Spectral Element Method
=================================================================

This directory contains the complete, organized MATLAB solution to 
1D spectral element method wave propagation.

DIRECTORY STRUCTURE:
--------------------

SEM1D/
├── wave.m                    : Main simulation program with real-time visualization
├── README.txt                : This file
├── output/                   : Directory for snapshot files
└── functions/                : All required library functions (12 files)
    ├── define_derivative_matrix.m
    ├── lagrange_deriv_GLL.m
    ├── zwgljd.m, zwgjd.m
    ├── jacg.m, jacobf.m
    ├── pnleg.m, pndleg.m
    ├── pnormj.m
    ├── gammaf.m
    ├── endw1.m, endw2.m

Note: The functions/ directory is automatically added to MATLAB's path when
wave.m is executed. No manual path management is required.

USAGE:
------

1. Running the simulation with real-time visualization:

   To run with Dirichlet boundary conditions (default):
   >> wave()
   or
   >> wave('dirichlet')

   To run with Neumann boundary conditions:
   >> wave('neumann')

SOLUTION IMPLEMENTATION:
------------------------

1. Global Mass Matrix (lines 117-124):
   - Assembles global mass matrix from local element contributions
   - Uses GLL quadrature weights, density, and Jacobian

2. Time Step Setup (lines 128-130):
   - Uses reference time step DT_REF = 0.25
   - Precomputes halfdt and halfdt2 for efficiency

3. Initial Conditions (lines 136-150):
   - Gaussian pulse centered at x = 50
   - Zero initial velocity and acceleration
   - Applies boundary conditions based on BC type

4. Predictor Step (lines 161-163):
   - Newmark time integration scheme
   - Updates displacement and velocity
   - Resets acceleration

5. Boundary Conditions (lines 167-174):
   - Dirichlet: fixes displacement at both ends to zero
   - Neumann: traction-free (natural) boundary conditions

6. Solver (lines 178-196):
   - Assembles global force vector from stiffness matrix
   - Stiffness computed using GLL quadrature and derivative matrix
   - Solves for acceleration using F = M * a

7. Corrector Step (lines 200-202):
   - Completes Newmark time integration
   - Updates acceleration and velocity
   - Displacement remains unchanged (already updated in predictor)

PARAMETERS:
-----------

Defined in wave.m:
- NSPEC = 20    : Number of spectral elements
- NGLL = 7      : Number of GLL points per element
- NGLOB = 121   : Total number of global points
- NSTEP = 2000  : Number of time steps
- DT = 0.25     : Time step in seconds
- LENGTH = 100  : Domain length
- DENSITY = 1   : Material density
- SHEARMODULUS = 1 : Shear modulus
