import numpy as np

# coefficients for the pressure equation at branching points, can choose gamma_2 and gamma_3 to 0 or 2
gamma_1 = 0 
gamma_2 = 2 
gamma_3 = 2

mu = 0.045 # dynamic viscosity 
rho = 1.050 # density

gamma_profile = 2 # coefficient of the a priory axial velocity profile
alpha = (gamma_profile + 2) / (gamma_profile + 1)
k_r = 2 * np.pi * (gamma_profile + 2) * mu / rho # friction factor