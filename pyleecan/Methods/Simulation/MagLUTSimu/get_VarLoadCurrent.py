from numpy import zeros

from ....Classes.VarLoadCurrent import VarLoadCurrent


def get_VarLoadCurrent(I_dq, speed):
    kwargs = {
        "name": "Variable Current",  # Name of the multi-simulation
        "desc": "Variable Current",  # Multi-simulation description
        "datakeeper_list": [],  # List containing DataKeepers to extract results
        "is_keep_all_output": True,  # store every output in a list
        "stop_if_error": True,  # Stop if a simulation fails
        # "ref_simu_index": 0,  # Index of the reference simulation, None if not part
        "nb_simu": 0,  # Number of simulations
        "is_reuse_femm_file": True,  # True to reuse the FEMM file for each simulation
        "postproc_list": [],  # post-proc. to run on XOutput after the multisimulation
        "pre_keeper_postproc_list": [],  # If not None, replace the reference simulation postproc_list in each generated simulation (run before datakeeper)
        "post_keeper_postproc_list": [],  # List of post-processing to run on output after each simulation (except reference one) after the datakeeper
        "type_OP_matrix": 1,  # Select with kind of OP_matrix is used 0: (N0,I0,Phi0,T,P), 1:(N0,Id,Iq,T,P)
        "is_torque": False,  # True if the Torque is defined in OP_matrix
        "is_power": False,  # True if the Power is defined in OP_matrix
    }
    var_load = VarLoadCurrent(**kwargs)

    # creating the Operating point matrix
    OP_matrix = zeros((I_dq.shape[0], 3))

    # set speed and current
    OP_matrix[:, 0] = speed
    OP_matrix[:, 1] = I_dq[:, 0]
    OP_matrix[:, 2] = I_dq[:, 1]

    # set variable parameter simulation
    var_load.OP_matrix = OP_matrix

    # return
    return var_load
