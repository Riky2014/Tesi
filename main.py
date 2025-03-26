from fenics import *

import time
import pickle
import logging
import numpy as np
from tqdm import tqdm

from geometry import branches, bifurcations
from solve import update_bcs, update_branching_points, solve_branches
from initialization import initialize_branches, create_mesh_and_function_spaces
from utils import plot_branches_combined, check_boundary_conditions, check_branching_consistency, check_inflow_data

# Set the logging level to suppress FFC messages
set_log_level(LogLevel.ERROR)
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL_LEGACY').setLevel(logging.WARNING)

# Check on network dictionaries
check_boundary_conditions(branches)
check_branching_consistency(branches, bifurcations)

# Simulation parameters
T = 5
T_base = 1

mult_fact = 8
h = mult_fact * 1 / 32
dt = mult_fact * 1e-5

N_T_base = int(round(T_base / dt))
N_T = int(round(T / dt)) + 1

print(f'\nh = {h} \ndt = {dt} \nCFL ratio : {h / dt}')
print()

# Define mesh, function space and initialize branches
branch_boundary_conditions = {}
create_mesh_and_function_spaces(branches, h)
initialize_branches(branches, branch_boundary_conditions)

# Define dict to store results
results_total = {branch_id: {'A': [], 'q': []} for branch_id in branches.keys()}

# Loop over all the time instant until final time T
for n in tqdm(range(N_T)):
     
    solve_branches(branches, branch_boundary_conditions, results_total, dt) 
    update_bcs(branches, branch_boundary_conditions, (n + 1) * dt, dt, mult_fact * (n % N_T_base)) 
    update_branching_points(bifurcations, branches, branch_boundary_conditions, dt)

    # if n % 100 == 0:
    #     plot_branches_combined(branches, n)

with open('closed_network_gamma2.pkl', 'wb') as f:
    pickle.dump(results_total, f)