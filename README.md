# Computational investigation of Circle of Willis hemodynamics by a one-dimensional model

## Overview
This project simulates blood flow through a vascular network using the finite element method (FEM) implemented in FEniCS. The model accounts for branching vessels and includes both physiological and pathological configurations. The simulation incorporates pressure dynamics at the branching points via a linearized Newton method and uses either area or flux inflow data.

## File Structure

- **Geometry**
    Contains the `.stl` file of the original 3D geometry and the `.vtp` of the 1D geometry, from where we can extract the lenght and reference area of each vessel. The `branching_information.ipynb` is used to compute the angles of each bifurcation

- **`geometry.py`**  
  Defines the vascular network geometry. It sets up:
  - **Branches:** Vessel properties (length, reference area, elasticity parameter) and boundary conditions.
  - **Bifurcations:** Details of network branch junctions and associated angles.
  - **Inflow Profiles:** Loads inflow data from CSV files for different inlet vessels.

- **`initialization.py`**  
  Sets up the simulation by:
  - Creating meshes and function spaces for each vessel using FEniCS.
  - Initializing flow variables (area and flux) on each branch based on boundary conditions and inflow data.
  - Establishing initial boundary condition dictionaries.

- **`parameters.py`**  
  Contains model parameters including:
  - Coefficients for the pressure equations at branching points (`gamma_1`, `gamma_2`, `gamma_3`).
  - Fluid properties such as dynamic viscosity (`mu`) and density (`rho`).
  - Parameters governing the axial velocity profile and friction factor.

- **`solve.py`**  
  Implements the core numerical routines:
  - **Finite Element Assembly:** Constructs matrices and vectors for the discretized governing equations.
  - **Local Solver:** Solves a linear variational problem on each branch.
  - **Newton Iteration at Branching Points:** Assembles and solves the linearized system to update pressure and flow at junctions.
  - **Boundary Condition Updates:** Computes and updates inlet/outlet boundary conditions for each time step.

- **`main.py`**  
  Acts as the driver script:
  - Validates network configuration and consistency using utility functions.
  - Initializes simulation parameters (time step, spatial resolution, total simulation time).
  - Calls routines to create meshes, initialize branches, and run the time-stepping loop.
  - Saves simulation results (areas and fluxes for each branch) into a pickle file for post-processing.

- **`utils.py`**  
  Provides helper functions to:
  - Check the validity of boundary conditions and the consistency of bifurcations.
  - Validate the inflow data type.
  - Plot area and flux profiles for branches to visualize simulation results.



