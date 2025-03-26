from fenics import *
from geometry import known_inflow_data, inflow_data_BAS, inflow_data_L_ICA, inflow_data_R_ICA

# Function to create meshes and function spaces
def create_mesh_and_function_spaces(branches, h):
    for branch_name, branch_data in branches.items():
        L = branch_data['L']
        N = int(L / h)
        mesh = IntervalMesh(N, 0, L)
        branch_data['mesh'] = mesh

        P1 = FiniteElement('P', mesh.ufl_cell(), 1)
        element = MixedElement([P1, P1])
        V = FunctionSpace(mesh, element)
        branch_data['V'] = V

        D1 = FiniteElement('DG', mesh.ufl_cell(), 0)
        element_der = MixedElement([D1, D1])
        V_der = FunctionSpace(mesh, element_der)
        branch_data['V_der'] = V_der


# Function to initialize branches and boundary conditions
def initialize_branches(branches, branch_boundary_conditions):
  for branch_id, branch_data in branches.items():
        A0 = branch_data['A0']

        if known_inflow_data == 'area':
          q0 = 0

        elif known_inflow_data == 'flux':

          if branch_id == '5':
            inflow_data = - inflow_data_L_ICA
          elif branch_id == '12':
            inflow_data = - inflow_data_R_ICA
          else:
            inflow_data = inflow_data_BAS

          q0 = inflow_data[0][1] * A0

        V = branch_data['V']

        uh_old = interpolate(Expression(('A0', 'q0'), degree=1, A0=A0, q0=q0), V)
        Ah_old, qh_old = uh_old.split()

        branch_data['uh'] = uh_old
        branch_data['Ah_old'] = Ah_old
        branch_data['qh_old'] = qh_old

        branch_boundary_conditions[branch_id] = {
        'left': {
            'A_left': A0,
            'q_left': q0
        },
        'right': {
            'A_right': A0,
            'q_right': q0
        }
       }