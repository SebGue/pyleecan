# -*- coding: utf-8 -*-
from os.path import join

from ....Classes.ElmerResultsVTU import ElmerResultsVTU


def get_meshsolution(self, output):
    """Build the MeshSolution objects from the FEA outputs.

    Parameters
    ----------
    self : MagElmer
        a MagElmer object
    output: Output
        An Output object

    Returns
    -------
    meshsol: MeshSolution
        a MeshSolution object with Elmer outputs at every time step
    """
    # logger
    logger = self.get_logger()

    # if meshsoltion is not requested set meshsolution to None
    if not self.is_get_mesh or not self.is_save_FEA:
        logger.info("StructElmer: MeshSolution is not stored by request.")
        output.struct.meshsolution = None
        return False

    # setup Elmer result helper class
    ElmerVtu = ElmerResultsVTU()

    ElmerVtu.store_dict = {
        "magnetic flux density e": {
            "name": "Magnetic Flux Density B",
            "unit": "T",
            "symbol": "B",
            "norm": 1,
        },
        "magnetic vector potential e": {
            "name": "Magnetic Vector Potential A",
            "unit": "Wb",
            "symbol": "A",
            "norm": 1,
        },
        "magnetic field strength e": {
            "name": "Magnetic Field H",
            "unit": "A/m",
            "symbol": "H",
            "norm": 1,
        },
        "current density e": {
            "name": "Current Density J",
            "unit": "A/mm2",
            "symbol": "J",
            "norm": 1,
        },
    }

    # ElmerVtu.store_dict = {
    #     "magnetic vector potential": {
    #         "name": "Magnetic Vector Potential A",
    #         "unit": "Wb",
    #         "symbol": "A",
    #         "norm": 1,
    #     },
    #     "magnetic flux density": {
    #         "name": "Magnetic Flux Density B",
    #         "unit": "T",
    #         "symbol": "B",
    #         "norm": 1,
    #     },
    #     "magnetic field strength": {
    #         "name": "Magnetic Field H",
    #         "unit": "A/m",
    #         "symbol": "H",
    #         "norm": 1,
    #     },
    #     "current density": {
    #         "name": "Current Density J",
    #         "unit": "A/mm2",
    #         "symbol": "J",
    #         "norm": 1,
    #     }
    # }

    ElmerVtu.label = "Elmer MagnetoDynamics"
    ElmerVtu.file_path = join(self.get_path_save_fea(output), "step_t0002.vtu")

    output.mag.meshsolution = ElmerVtu.build_meshsolution(is_point_data=False)

    return True
