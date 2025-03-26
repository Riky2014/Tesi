import numpy as np
from fenics import *
from parameters import gamma_1, gamma_2, gamma_3, rho, alpha, k_r
from geometry import known_inflow_data, inflow_data_BAS, inflow_data_L_ICA, inflow_data_R_ICA


# Assemble matrix and vectors for the FE discretization
def H(A, q, A0, beta):
  return as_tensor([[0, 1], [beta / (2 * rho * A0) * A ** 0.5 - alpha * (q / A) ** 2, 2 * alpha * q / A]])

def F(A, q, A0, beta):
  return as_vector([q, beta / (3 * rho * A0) * A ** 1.5 - beta / (3 * rho) * A0 ** 0.5 + alpha * q ** 2 / A])

def B(A, q):
  return as_vector([0, k_r * q / A])

def S(A, q):
  return as_vector([0, k_r * q / A])

def dS_dU(A, q):
  return as_tensor([[0, 0], [- k_r * q / A ** 2, k_r / A]])

# Define where the left and right boundary are
def boundary_L(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol)

def boundary_R(x, on_boundary):
    tol = 1E-14
    return on_boundary and x[0]>0.1

# Solve the variational problem for a single branch 
def LinearProblem(bc, V, A_old, q_old, A0, beta, dt):
  # Define trial functions and test functions
  A, q = TrialFunctions(V)
  v, z = TestFunctions(V)
  uh = Function(V)

  a = inner(A, v) * dx + inner(q, z) * dx

  l = (
        A_old * v * dx
      + q_old * z * dx
      + dt * ((F(A_old, q_old, A0, beta) - dt / 2 * dot(H(A_old, q_old, A0, beta), S(A_old, q_old))))[0] * v.dx(0) * dx
      + dt * ((F(A_old, q_old, A0, beta) - dt / 2 * dot(H(A_old, q_old, A0, beta), S(A_old, q_old))))[1] * z.dx(0) * dx
      + dt ** 2 / 2 * (dot(dS_dU(A_old, q_old), F(A_old, q_old, A0, beta).dx(0)))[0] * v * dx
      + dt ** 2 / 2 * (dot(dS_dU(A_old, q_old), F(A_old, q_old, A0, beta).dx(0)))[1] * z * dx
      - dt ** 2 / 2 * (dot(H(A_old, q_old, A0, beta), F(A_old, q_old, A0, beta).dx(0)))[0] * v.dx(0) * dx
      - dt ** 2 / 2 * (dot(H(A_old, q_old, A0, beta), F(A_old, q_old, A0, beta).dx(0)))[1] * z.dx(0) * dx
      - dt * (S(A_old, q_old) - dt / 2 * dot(dS_dU(A_old, q_old), S(A_old, q_old)))[0] * v * dx
      - dt * (S(A_old, q_old) - dt / 2 * dot(dS_dU(A_old, q_old), S(A_old, q_old)))[1] * z * dx
  )

  solve(a == l, uh, bc)

  return uh

# Loop over all the branches to find the new solution and update the old one
def solve_branches(branches, branch_boundary_conditions, results_total, dt):
    for branch_id, data in branches.items():
        # Extract boundary conditions
        bc_data = branch_boundary_conditions[branch_id]

        A_left = bc_data['left']['A_left']
        q_left = bc_data['left']['q_left']
        A_right = bc_data['right']['A_right']
        q_right = bc_data['right']['q_right']

        # Initialize boundary conditions
        bcs = []

        bc_A_left = DirichletBC(data['V'].sub(0), A_left, boundary_L)
        bc_q_left = DirichletBC(data['V'].sub(1), q_left, boundary_L)
        bc_A_right = DirichletBC(data['V'].sub(0), A_right, boundary_R)
        bc_q_right = DirichletBC(data['V'].sub(1), q_right, boundary_R)
        
        bcs.extend([bc_A_left, bc_q_left, bc_A_right, bc_q_right])

        # Solve the local problem
        uh = LinearProblem(bcs, data['V'], data['Ah_old'], data['qh_old'], data['A0'], data['beta'], dt)
        Ah, qh = uh.split(deepcopy=True)
        data['uh'] = uh  # Store the computed solution
        data['Ah_old'].assign(Ah)
        data['qh_old'].assign(qh)

        # Store the results
        results_total[branch_id]['A'].append(Ah.compute_vertex_values(data['mesh']))
        results_total[branch_id]['q'].append(qh.compute_vertex_values(data['mesh']))


# Those are necessary to compute pointwise values for boundary condition
def U(A, q):
  return np.array([A, q])

def H_vec(A, q, A0, beta):
  return np.array([[0, 1], [beta / (2 * rho * A0) * A ** 0.5 - alpha * (q / A) ** 2, 2 * alpha * q / A]])

def B_vec(A, q):
  return np.array([0, k_r * q / A])

def S_vec(A, q):
  return np.array([0, k_r * q / A])

def c_alpha(A, q, A0, beta):
  return (beta / (2 * rho * A0) * A ** 0.5 + (q / A) ** 2 * alpha * (alpha - 1)) ** 0.5

def l1(A, q, A0, beta):
  return np.array([c_alpha(A, q, A0, beta) - alpha * q / A, 1.])

def l2(A, q, A0, beta):
  return np.array([- c_alpha(A, q, A0, beta) - alpha * q / A, 1.])

def CC(A, q, A0, beta, u_der, V_der, x, dt):
  du_dz = project(u_der, V_der)(x)
  return U(A, q) - dt * H_vec(A, q, A0, beta) @ du_dz - dt * B_vec(A, q)


# Impose Dirichlet boundary condition at the inlet
def inlet_bc_q(A, q, A0, beta, u_der, V_der, x, dt):
  if x < 1e-6:
    q_inlet = (np.dot(l2(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l2(A(x), q(x), A0, beta)[0] * U(A(x), q(x))[0] ) / l2(A(x), q(x), A0, beta)[1]
  else:
    q_inlet = (np.dot(l1(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l1(A(x), q(x), A0, beta)[0] * U(A(x), q(x))[0] ) / l1(A(x), q(x), A0, beta)[1]

  return q_inlet

def inlet_bc_A(A, q, A0, beta, u_der, V_der, x, dt):
  if x < 1e-6:
    A_inlet = (np.dot(l2(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l2(A(x), q(x), A0, beta)[1] * U(A(x), q(x))[1] ) / l2(A(x), q(x), A0, beta)[0]
  else:
     A_inlet = (np.dot(l1(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l1(A(x), q(x), A0, beta)[1] * U(A(x), q(x))[1] ) / l1(A(x), q(x), A0, beta)[0]

  return A_inlet

# Impose constant P0 at the outlet
def outlet_bc(A, q, A0, beta, u_der, V_der, x, dt):
  if x < 1e-6:
    q_inlet = (np.dot(l2(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l2(A(x), q(x), A0, beta)[0] * U(A(x), q(x))[0] ) / l2(A(x), q(x), A0, beta)[1]
  else:
    q_inlet = (np.dot(l1(A(x), q(x), A0, beta), CC(A(x), q(x), A0, beta, u_der, V_der, x, dt)) - l1(A(x), q(x), A0, beta)[0] * U(A(x), q(x))[0] ) / l1(A(x), q(x), A0, beta)[1]

  return A0, q_inlet


# Update inlet and outlet boundary constion for the new time step
def update_bcs(branches, branch_boundary_conditions, t, dt, n):
    for branch_id, data in branches.items():

      # Retrieve the boundary conditions for the start of the branch
      left_boundary = data['boundary']['left']

      if left_boundary['type'] == 'inlet':

        if known_inflow_data == 'area':
          A_inlet = data['A0'] + 0.01 * t
          q_inlet = inlet_bc_q(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], 0, dt)

        elif known_inflow_data == 'flux':

          if branch_id == '1':
            inflow_data = inflow_data_BAS
          elif branch_id == '5':
            inflow_data = inflow_data_L_ICA
          elif branch_id == '12':
            inflow_data = inflow_data_R_ICA

          q_inlet = inflow_data[n][1] * branch_boundary_conditions[branch_id]['left']['A_left']
          A_inlet = inlet_bc_A(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], 0, dt)

        branch_boundary_conditions[branch_id]['left'] = {'A_left': A_inlet, 'q_left': q_inlet}

      elif left_boundary['type'] == 'outlet':
        A_outlet, q_outlet = outlet_bc(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], 0, dt)
        branch_boundary_conditions[branch_id]['left'] = {'A_left': A_outlet, 'q_left': q_outlet}

      # Retrieve the boundary conditions for the end of the branch
      right_boundary = data['boundary']['right']

      if right_boundary['type'] == 'inlet':

        if known_inflow_data == 'area':
          A_inlet = data['A0'] + 0.01 * t
          q_inlet = inlet_bc_q(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], data['L'], dt)

        elif known_inflow_data == 'flux':

          if branch_id == '1':
            inflow_data = inflow_data_BAS
          elif branch_id == '5':
            inflow_data = inflow_data_L_ICA
          elif branch_id == '12':
            inflow_data = inflow_data_R_ICA

          q_inlet = - inflow_data[n][1] * branch_boundary_conditions[branch_id]['right']['A_right']
          A_inlet = inlet_bc_A(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'],  data['L'], dt)

        branch_boundary_conditions[branch_id]['right'] = {'A_right': A_inlet, 'q_right': q_inlet}

      elif right_boundary['type'] == 'outlet':
        A_outlet, q_outlet = outlet_bc(data['Ah_old'], data['qh_old'], data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], data['L'], dt)
        branch_boundary_conditions[branch_id]['right'] = {'A_right': A_outlet, 'q_right': q_outlet}


# Define terms of the linearized newton matrix
def C_A_in(A, q, A0, beta, gamma):
  return beta / (2 * rho * A0 * sqrt(A)) - (q ** 2 / A ** 3) * (1 - 2 * gamma * np.sign(q))

def C_q_in(A, q, A0, gamma):
  return (q / A ** 2) * (1 - 2 * gamma * np.sign(q))

def b_in(A, q, A0, beta, gamma):
  return beta / rho * (A - 2 * sqrt(A * A0)) / (2 * A0 * sqrt(A)) + 0.5 * (q / A) ** 2 * (1 - 2 * gamma * np.sign(q))

def C_A_out(A, q, A0, theta, beta, gamma):
  return beta / (2 * rho * A0 * sqrt(A)) - (q ** 2 / A ** 3) * (1 + 2 * gamma * sqrt(2 * (1 - cos(theta))) * np.sign(q))

def C_q_out(A, q, A0, theta, gamma):
  return (q / A ** 2) * (1 + 2 * gamma * sqrt(2 * (1 - cos(theta))) * np.sign(q))

def b_out(A, q, A0, theta, beta, gamma):
  return beta / rho * (A - 2 * sqrt(A * A0)) / (2 * A0 * sqrt(A)) + 0.5 * (q / A) ** 2 * (1 + 2 * gamma * sqrt(2 * (1 - cos(theta))) * np.sign(q))


# Check if the linearized solution solve 2 non linear equations
def g_1(A, q, gamma):
  return gamma * (q / A) ** 2

def g_i(A, q, gamma, theta):
  return gamma * (q / A) ** 2 * sqrt(2 * (1 - cos(theta)))

def check_non_linear_system(A1, A2, A3, q1, q2, q3, data, connected_data1, connected_data2, bifurcation):
    theta_2 = bifurcation['theta'][1]
    theta_3 = bifurcation['theta'][2]

    res_2 = data['beta'] / (rho * data['A0']) * (np.sqrt(A1) - np.sqrt(data['A0'])) + 0.5 * (q1 / A1) ** 2 - np.sign(q1 / A1) * g_1(A1, q1, gamma_1) \
            - connected_data1['beta'] / (rho * connected_data1['A0']) * (np.sqrt(A2) - np.sqrt(connected_data1['A0'])) - 0.5 * (q2 / A2) ** 2 - np.sign(q2 / A2) * g_i(A2, q2, gamma_2, theta_2)

    res_3 = data['beta'] / (rho * data['A0']) * (np.sqrt(A1) - np.sqrt(data['A0'])) + 0.5 * (q1 / A1) ** 2 - np.sign(q1 / A1) * g_1(A1, q1, gamma_1) \
            - connected_data2['beta'] / (rho * connected_data2['A0']) * (np.sqrt(A3) - np.sqrt(connected_data2['A0'])) - 0.5 * (q3 / A3) ** 2 - np.sign(q3 / A3) * g_i(A3, q3, gamma_3, theta_3)

    return res_2, res_3

# Assemble and solve the linearized newton system for each branch
def assemble_solve_newton(data, connected_data1, connected_data2, bifurcation, positions, dt, toll=1e-5):
    # Define positions
    x1 = data['L'] if positions[0] == 'right' else 0
    x2 = connected_data1['L'] if positions[1] == 'right' else 0
    x3 = connected_data2['L'] if positions[2] == 'right' else 0

    # Get old values
    A_branch_old_1 = data['Ah_old'](x1)
    A_branch_old_2 = connected_data1['Ah_old'](x2)
    A_branch_old_3 = connected_data2['Ah_old'](x3)

    q_branch_old_1 = data['qh_old'](x1)
    q_branch_old_2 = connected_data1['qh_old'](x2)
    q_branch_old_3 = connected_data2['qh_old'](x3)

    # Initialize old values for new branches
    A_branch_0_1, A_branch_0_2, A_branch_0_3 = A_branch_old_1, A_branch_old_2, A_branch_old_3
    q_branch_0_1, q_branch_0_2, q_branch_0_3 = q_branch_old_1, q_branch_old_2, q_branch_old_3

    # theta values
    theta_2 = bifurcation['theta'][1]
    theta_3 = bifurcation['theta'][2]

    err = float('inf')
    iter = 0

    while (err > toll) and (iter < 100):
        iter += 1

        lhs = np.zeros((6, 6))
        rhs = np.zeros(6)

        if positions[0] == 'right':
          lhs[0, 3] = 1
        else: 
          lhs[0, 3] = -1
        if positions[1] == 'right':
          lhs[0, 4] = 1
        else: 
          lhs[0, 4] = -1
        if positions[2] == 'right':
          lhs[0, 5] = 1
        else: 
          lhs[0, 5] = -1  

        lhs[1, 0] = C_A_in(A_branch_old_1, q_branch_old_1, data['A0'], data['beta'], gamma_1)
        lhs[1, 1] = -C_A_out(A_branch_old_2, q_branch_old_2, connected_data1['A0'], theta_2, connected_data1['beta'], gamma_2)
        lhs[1, 3] = C_q_in(A_branch_old_1, q_branch_old_1, data['A0'], gamma_1)
        lhs[1, 4] = -C_q_out(A_branch_old_2, q_branch_old_2, connected_data1['A0'], theta_2, gamma_2)

        lhs[2, 0] = C_A_in(A_branch_old_1, q_branch_old_1, data['A0'], data['beta'], gamma_1)
        lhs[2, 2] = -C_A_out(A_branch_old_3, q_branch_old_3, connected_data2['A0'], theta_3, connected_data2['beta'], gamma_3)
        lhs[2, 3] = C_q_in(A_branch_old_1, q_branch_old_1, data['A0'], gamma_1)
        lhs[2, 5] = -C_q_out(A_branch_old_3, q_branch_old_3, connected_data2['A0'], theta_3, gamma_3)

        if positions[0] == 'right':
          lhs[3, 0] = l1(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'])[0]
          lhs[3, 3] = l1(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'])[1]
        else:
          lhs[3, 0] = l2(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'])[0]
          lhs[3, 3] = l2(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'])[1]

        if positions[1] == 'right':
          lhs[4, 1] = l1(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'])[0]
          lhs[4, 4] = l1(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'])[1]
        else:
          lhs[4, 1] = l2(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'])[0]
          lhs[4, 4] = l2(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'])[1]

        if positions[2] == 'right':
          lhs[5, 2] = l1(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'])[0]
          lhs[5, 5] = l1(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'])[1]
        else:
          lhs[5, 2] = l2(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'])[0]
          lhs[5, 5] = l2(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'])[1]

        rhs[1] = b_out(A_branch_old_2, q_branch_old_2, connected_data1['A0'], theta_2, connected_data1['beta'], gamma_2) - b_in(A_branch_old_1, q_branch_old_1, data['A0'], data['beta'], gamma_1)
        rhs[2] = b_out(A_branch_old_3, q_branch_old_3, connected_data2['A0'], theta_3, connected_data2['beta'], gamma_3) - b_in(A_branch_old_1, q_branch_old_1, data['A0'], data['beta'], gamma_1)

        if positions[0] == 'right':
            rhs[3] = np.dot(l1(A_branch_0_1, q_branch_0_1, data['A0'], data['beta']), CC(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], x1, dt))
        else:
            rhs[3] = np.dot(l2(A_branch_0_1, q_branch_0_1, data['A0'], data['beta']), CC(A_branch_0_1, q_branch_0_1, data['A0'], data['beta'], data['uh'].dx(0), data['V_der'], x1, dt))

        if positions[1] == 'right':
            rhs[4] = np.dot(l1(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta']), CC(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'], connected_data1['uh'].dx(0), connected_data1['V_der'], x2, dt))
        else:
            rhs[4] = np.dot(l2(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta']), CC(A_branch_0_2, q_branch_0_2, connected_data1['A0'], connected_data1['beta'], connected_data1['uh'].dx(0), connected_data1['V_der'], x2, dt))

        if positions[2] == 'right':
            rhs[5] = np.dot(l1(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta']), CC(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'], connected_data2['uh'].dx(0), connected_data2['V_der'], x3, dt))
        else:
            rhs[5] = np.dot(l2(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta']), CC(A_branch_0_3, q_branch_0_3, connected_data2['A0'], connected_data2['beta'], connected_data2['uh'].dx(0), connected_data2['V_der'], x3, dt))


        A_branch_1, A_branch_2, A_branch_3, q_branch_1, q_branch_2, q_branch_3 = np.linalg.solve(lhs, rhs)

        err = abs(A_branch_1 - A_branch_old_1) / A_branch_0_1 + abs(A_branch_2 - A_branch_old_2) / A_branch_0_2 + abs(A_branch_3 - A_branch_old_3) / A_branch_0_3 \
            + abs(q_branch_1 - q_branch_old_1) / q_branch_0_1 + abs(q_branch_2 - q_branch_old_2) / q_branch_0_2 + abs(q_branch_3 - q_branch_old_3) / q_branch_0_3

        A_branch_old_1, A_branch_old_2, A_branch_old_3 = A_branch_1, A_branch_2, A_branch_3
        q_branch_old_1, q_branch_old_2, q_branch_old_3 = q_branch_1, q_branch_2, q_branch_3

    # res_2, res_3 = check_non_linear_system(
    #     A_branch_1, A_branch_2, A_branch_3,
    #     q_branch_1, q_branch_2, q_branch_3,
    #     data, connected_data1, connected_data2, bifurcation
    # )

    #print(res_2, res_3)

    return A_branch_1, A_branch_2, A_branch_3, q_branch_1, q_branch_2, q_branch_3, err, iter


# Loop all over the bifurcations to update A and q  at branching points
def update_branching_points(bifurcations, branches, branch_boundary_conditions, dt):
    for bif_id, bif_data in bifurcations.items():
        branch_id = bif_data['branch_id'][0]  # Main branch
        connected_branches = bif_data['branch_id'][1:]  # Connected branches

        if len(connected_branches) == 2:
            cb1_id, cb2_id = connected_branches
            main_data = branches[branch_id]
            connected_data1 = branches[cb1_id]
            connected_data2 = branches[cb2_id]
            positions = bif_data['positions']

            A_branch_1, A_branch_2, A_branch_3, q_branch_1, q_branch_2, q_branch_3, err, iter = assemble_solve_newton(
                main_data, connected_data1, connected_data2, bif_data, positions, dt, toll=1e-5
            )

            # Update boundary conditions
            branches_to_update = [branch_id, cb1_id, cb2_id]
            values = [(A_branch_1, q_branch_1), (A_branch_2, q_branch_2), (A_branch_3, q_branch_3)]

            for branch, (A, q), pos in zip(branches_to_update, values, positions):
                branch_boundary_conditions[branch][pos] = {'A_' + pos: A, 'q_' + pos: q}

        else:
            raise ValueError("A branching point should connect to exactly two branches.")