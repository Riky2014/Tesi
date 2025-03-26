import numpy as np

# CLOSED network geometric configuration (Physiological conditions)

# Branch dictionary
# L: vessel lenght
# A0: reference area
# beta: elasticity parameter
# boundary: type of boundary condition at left and right hand

branches = {
    '1': {
        'L': 1.28,
        'A0': 0.0740,
        'beta': 0.0725e7,
        'boundary': {
            'left': {'type': 'inlet'},
            'right': {'type': 'branch'}
        }
    },
    '2': {
        'L': 1.3212,
        'A0': 0.0347,
        'beta': 0.0993e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '3': {
        'L': 3.0909,
        'A0': 0.0338,
        'beta': 0.0981e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
    '4': {
        'L': 1.1713,
        'A0': 0.0256,
        'beta': 0.0853e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '5': {
        'L': 4.2887,
        'A0': 0.1252,
        'beta': 0.0944e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'inlet'}
        }
    },
    '6': {
        'L': 0.4203,
        'A0': 0.0756,
        'beta': 0.1467e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '7': {
        'L': 3.0665,
        'A0': 0.0288,
        'beta': 0.0905e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '8': {
        'L': 4.4854,
        'A0': 0.0580,
        'beta': 0.1284e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
    '9': {
        'L': 1.0012,
        'A0': 0.0295,
        'beta': 0.0916e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '10': {
        'L': 2.5496,
        'A0': 0.0331,
        'beta': 0.0970e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
    '11': {
        'L': 1.8936,
        'A0': 0.0301,
        'beta': 0.0926e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '12': {
        'L': 4.8650,
        'A0': 0.1260,
        'beta': 0.0946e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'inlet'}
        }
    },
    '13': {
        'L': 0.6961,
        'A0': 0.0800,
        'beta': 0.1509e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '14': {
        'L': 2.1127,
        'A0': 0.0435,
        'beta': 0.1113e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '15': {
        'L': 3.3325,
        'A0': 0.0623,
        'beta': 0.1331e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
    '16': {
        'L': 0.4632,
        'A0': 0.0214,
        'beta': 0.0781e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'branch'}
        }
    },
    '17': {
        'L': 4.4019,
        'A0': 0.0281,
        'beta': 0.0894e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
    '18': {
        'L': 3.9169,
        'A0': 0.0323,
        'beta': 0.0958e7,
        'boundary': {
            'left': {'type': 'branch'},
            'right': {'type': 'outlet'}
        }
    },
}

# Bifurcations dict
# branch_id: ids of the branches involved in the bifurcation
# positions: where the bif is located for each branch
# theta: angle between the first branch and each of the other 2

bifurcations = {
    'bif_1': {
        'branch_id': ['1', '2', '9'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 0.7748, 0.7771]
    },
    'bif_2': {
        'branch_id': ['2', '3', '4'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 0.7054, 0.9872]
    },
    'bif_3': {
        'branch_id': ['5', '4', '6'],
        'positions': ['left', 'right', 'left'],
        'theta': [None, 1.2182, 0.5778]
    },
    'bif_4': {
        'branch_id': ['6', '7', '8'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 1.9317, 0.1917]
    },
    'bif_5': {
        'branch_id': ['7', '16', '17'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 2.4837, 0.2554]
    },
    'bif_6': {
        'branch_id': ['9', '10', '11'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 0.8728, 1.1379]
    },
    'bif_7': {
        'branch_id': ['12', '13', '11'],
        'positions': ['left', 'left', 'right'],
        'theta': [None, 0.6948, 1.3204]
    },
    'bif_8': {
        'branch_id': ['13', '14', '15'],
        'positions': ['right', 'left', 'left'],
        'theta': [None, 1.0024, 0.9820]
    },
    'bif_9': {
        'branch_id': ['14', '16', '18'],
        'positions': ['right', 'right', 'left'],
        'theta': [None, 0.8910, 0.9184]
    }
}

# Inflow profiles for each of the inlet vessel
# inflow data: flux or area
# BAS: basilar artery
# L_ICA: Left Internal Carotid artery
# R_ICA: Right Internal Carotid artery

known_inflow_data = 'flux'
inflow_data_BAS = np.loadtxt('Inflow_profiles/BAS_inflow.csv', delimiter=',', skiprows=1)
inflow_data_L_ICA = np.loadtxt('Inflow_profiles/L_ICA_inflow.csv', delimiter=',', skiprows=1)
inflow_data_R_ICA = np.loadtxt('Inflow_profiles/R_ICA_inflow.csv', delimiter=',', skiprows=1)

# # TREE network geometric configuration (Patological conditions: vessel 7 is missing)

# branches = {
#     '1': { #0
#         'L': 1.28,
#         'A0': 0.0740,
#         'beta': 0.0725e7,
#         'boundary': {
#             'left': {'type': 'inlet'},
#             'right': {'type': 'branch'}
#         }
#     },
#     '2': {#2
#         'L': 1.3212,
#         'A0': 0.0347,
#         'beta': 0.0993e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '3': {#4
#         'L': 3.0909,
#         'A0': 0.0338,
#         'beta': 0.0981e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },

#     '4': {#5
#         'L': 1.1713,
#         'A0': 0.0256,
#         'beta': 0.0853e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '5': {#7
#         'L': 4.2887,
#         'A0': 0.1252,
#         'beta': 0.0944e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'inlet'}
#         }
#     },

#     '21': { #8 (6+8)
#         'L': 4.9727,
#         'A0': 0.0602,
#         'beta': 0.1309e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },

#     '9': {#9
#         'L': 1.0012,
#         'A0': 0.0295,
#         'beta': 0.0916e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '10': {#21
#         'L': 2.5496,
#         'A0': 0.0331,
#         'beta': 0.0970e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },

#     '11': {#11
#         'L': 1.8936,
#         'A0': 0.0301,
#         'beta': 0.0926e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '12': {#20
#         'L': 4.8650,
#         'A0': 0.1260,
#         'beta': 0.0946e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'inlet'}
#         }
#     },

#     '13': {#13
#         'L': 0.6961,
#         'A0': 0.0800,
#         'beta': 0.1509e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '14': {#15
#         'L': 2.1127,
#         'A0': 0.0435,
#         'beta': 0.1113e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'branch'}
#         }
#     },

#     '15': {#19
#         'L': 3.3325,
#         'A0': 0.0623,
#         'beta': 0.1331e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },

#     '22': { #17 (16+17)
#         'L': 4.1224,
#         'A0': 0.0271,
#         'beta': 0.0878e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },

#     '18': {#18
#         'L': 3.9169,
#         'A0': 0.0323,
#         'beta': 0.0958e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },
# }


# # Define the bifurcations dictionary
# bifurcations = {
#     'bif_1': {
#         'branch_id': ['1', '2', '9'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, 0.7748, 0.7771]
#     },
#     'bif_2': {
#         'branch_id': ['2', '3', '4'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, 0.7054, 0.9872]
#     },
#     'bif_3': {
#         'branch_id': ['5', '4', '21'],
#         'positions': ['left', 'right', 'left'],
#         'theta': [None, 1.2182, 0.5778]
#     },
#     'bif_4': {
#         'branch_id': ['9', '10', '11'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, 0.8728, 1.1379]
#     },
#     'bif_5': {
#         'branch_id': ['12', '13', '11'],
#         'positions': ['left', 'left', 'right'],
#         'theta': [None, 0.6948, 1.3204]
#     },
#     'bif_6': {
#         'branch_id': ['13', '14', '15'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, 1.0024, 0.9820]
#     },
#     'bif_7': {
#         'branch_id': ['14', '22', '18'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, 0.8910, 0.9184]
#     }
# }

# # Inflow profiles for each of the inlet vessel
# # inflow data: flux or area
# # BAS: basilar artery
# # L_ICA: Left Internal Carotid artery
# # R_ICA: Right Internal Carotid artery

# known_inflow_data = 'flux'
# inflow_data_BAS = np.loadtxt('Inflow_profiles/BAS_inflow.csv', delimiter=',', skiprows=1)
# inflow_data_L_ICA = np.loadtxt('Inflow_profiles/L_ICA_inflow.csv', delimiter=',', skiprows=1)
# inflow_data_R_ICA = np.loadtxt('Inflow_profiles/R_ICA_inflow.csv', delimiter=',', skiprows=1)


# # SINGLE BRANCH network used for convergence test

# known_inflow_data = 'area'

# branches = {
#     '1': {
#         'L': 1,
#         'A0': 0.126,
#         'beta': 0.060606e7,
#         'boundary': {
#             'left': {'type': 'inlet'},
#             'right': {'type': 'branch'}
#         }
#     },
#     '2': {
#         'L': 1,
#         'A0': 0.126,
#         'beta': 0.060606e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },
#     '3': {
#         'L': 2,
#         'A0': 0.126,
#         'beta': 0.060606e7,
#         'boundary': {
#             'left': {'type': 'branch'},
#             'right': {'type': 'outlet'}
#         }
#     },
# }


# # Define the bifurcations dictionary
# bifurcations = {
#     'bif_1': {
#         'branch_id': ['1', '2', '3'],
#         'positions': ['right', 'left', 'left'],
#         'theta': [None, np.pi / 3, np.pi / 4]
#     }
# }