import matplotlib.pyplot as plt

# checks if for every vessel we have imposed valid boundary conditions, that are {inlet, outlet, branch}
def check_boundary_conditions(branches):
    allowed_boundary_types = {'inlet', 'outlet', 'branch'}
    errors = []

    for branch_id, data in branches.items():
        for side, boundary in data['boundary'].items():
            if boundary['type'] not in allowed_boundary_types:
                errors.append(f"Error in branch '{branch_id}' on '{side}' side: Invalid boundary type '{boundary['type']}'.")

    if errors:
        for error in errors:
            print(error)
    else:
        print("All boundary conditions are valid.")


# checks if the branches and bifurcations dict are consistent
def check_branching_consistency(branches, bifurcations):
    # Create a dictionary to map branch IDs to their boundary types
    branch_boundary_types = {}
    for branch_id, data in branches.items():
        for side, boundary in data['boundary'].items():
            if boundary['type'] == 'branch':
                if branch_id not in branch_boundary_types:
                    branch_boundary_types[branch_id] = {}
                branch_boundary_types[branch_id][side] = True

    # Track whether any inconsistencies are found
    inconsistencies_found = False

    # Check the bifurcation dictionary for consistency
    for bif_id, bif_data in bifurcations.items():
        branch_ids = bif_data['branch_id']
        positions = bif_data['positions']

        for i, branch_id in enumerate(branch_ids):
            if branch_id in branch_boundary_types:
                side_positions = branch_boundary_types[branch_id]

                if positions[i] not in side_positions:
                    print(f"Inconsistency found in {bif_id}: Branch '{branch_id}' expected at '{list(side_positions.keys())}' but found at '{positions[i]}'.")
                    inconsistencies_found = True
                else:
                    # If matched, remove the checked side to avoid double-checking
                    del side_positions[positions[i]]

    # Print a message if all bifurcations are consistent
    if not inconsistencies_found:
        print("All bifurcations are consistent.")


# chevalidcks if the selected inflow data is 
def check_inflow_data(known_inflow_data):
    acceptable_inflow_data = ['area', 'flux']

    if known_inflow_data not in acceptable_inflow_data:
        print(f"Invalid known inflow data: '{known_inflow_data}'. It must be either 'area' or 'flux'.")
    else:
        print(f"Known inflow data '{known_inflow_data}' is valid.")


# plot area and flux of each vessel at a given time instant
def plot_branches_combined(branches, n, num_cols=5):
    num_branches = len(branches)
    num_rows = (num_branches + num_cols - 1) // num_cols  # Calculate required rows

    # Create a figure with subplots for both area and flux
    fig, axes = plt.subplots(num_rows * 2, num_cols, figsize=(num_cols * 3.8, num_rows * 8))
    axes = axes.flatten()

    for i, (branch_name, branch_data) in enumerate(branches.items()):
        if i >= num_cols * num_rows:  # Ensure we do not exceed the available axes
            break

        # Plot area in red on the first row
        ax_area = axes[i]
        Ah_old = branch_data['Ah_old']
        mesh = branch_data['mesh']
        ax_area.plot(mesh.coordinates(), Ah_old.compute_vertex_values(mesh), 'r', label='Area')
        ax_area.set_title(f'{branch_name}')
        ax_area.legend()

        # Plot flux in blue on the second row
        ax_flux = axes[num_cols * num_rows + i]
        qh_old = branch_data['qh_old']
        ax_flux.plot(mesh.coordinates(), qh_old.compute_vertex_values(mesh), 'b', label='Flux')
        ax_flux.legend()

    # Hide any unused subplots
    for j in range(i + 1, num_cols * num_rows):
        fig.delaxes(axes[j])
        fig.delaxes(axes[num_cols * num_rows + j])

    # Set global title
    fig.suptitle(f'Area and Flux Profile at Time Step {n}', fontsize=16)

    plt.tight_layout(rect=[0, 0, 1, 0.92])  # Adjust layout to make room for the global title
    plt.show()