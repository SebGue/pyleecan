from numpy import ones, pi, array, linspace, newaxis
from multiprocessing import cpu_count

from ....Classes.MagLUT import MagLUT
from ....Classes.InputCurrent import InputCurrent
from ....Classes.ImportMatrixVal import ImportMatrixVal
from ....Classes.MagFEMM import MagFEMM
from ....Classes.DataKeeper import DataKeeper


default_kwargs = {
    "type_BH_stator": 0,  # 0 == B(H) curve, 1 == mur_lin, 2 == inf. perm.
    "type_BH_rotor": 0,  # 0 == B(H) curve, 1 == mur_lin, 2 == inf. perm.
    "angle_stator_shift": 0,  # Angular position shift of the stator
    "file_name": "",  # Name of the file to save the FEMM model
    "Kmesh_fineness": 1,  # mesh fineness (1:default ,>1: finner ,<1: less fine)
    "Kgeo_fineness": 1,  # geometry fineness (1:default ,>1: finner ,<1: less fine)
    "is_periodicity_a": True,  # use machine symmetry
    "is_periodicity_t": False,  # use time symmetry
    "is_get_meshsolution": True,  # To get FEA mesh for latter post-procesing
    "is_save_meshsolution_as_file": False,  # To save FEA results in a dat file
    "nb_worker": cpu_count(),  # Number of threads for parallel FEMM execution
}

# TODO maybe do 'external' ref_simu so there is no need to re/store ref. machine
def run_simu(self):
    """Run the CurrentCharacteristic Simulation"""
    # store original machine to set it back after simulation
    org_machine = self.machine

    # normalize the machine
    self.machine = self.get_norm_machine(self.machine)

    # get some values for the simulation setup
    periodicity, is_anti_per = self.machine.comp_periodicity()[0:2]
    # periodicity = periodicity * (1 + is_anti_per)
    Nrev = self.Nrev
    Zs = self.machine.stator.slot.Zs
    Zs_per = Zs // periodicity

    Nr = Zs_per * 4 if self.Nr is None else self.Nr

    if Nrev is None:
        Nrev = 1 / periodicity  # rotate over one sym. section

    if Nrev == 0:
        Nr = 1

    # simulation settings
    kwargs = default_kwargs
    if self.simu_dict is not None:
        for key in self.simu_dict.keys():
            kwargs[key] = self.simu_dict[key]

    Na = 512 * self.machine.get_pole_pair_number()

    # setup the simulation output
    if self.parent is None or not isinstance(self.parent, MagLUT):
        output = MagLUT(simu=self)
    else:
        output = self.parent

    # set name
    if not self.name:
        self.name = self.machine.name + "MagLUT"

    # Defining magnetic simulation and input
    self.input = InputCurrent()
    self.mag = MagFEMM(**kwargs)

    # time discretization [s]
    # use number of revolutions and time steps instead of time object
    self.input.Nrev = Nrev
    self.input.Nt_tot = Nr

    # Angular discretization along the airgap circumference for flux density calculation
    angle = linspace(start=0, stop=2 * pi, num=Na, endpoint=False)
    self.input.angle = ImportMatrixVal()
    self.input.angle.value = angle

    # get current density
    S_ref = self.S_ref

    # compute the reference currents
    S_wind_act = self.machine.stator.winding.conductor.comp_surface_active()
    I_dq = array(S_ref * S_wind_act)

    if len(I_dq.shape) == 1:
        I_dq = I_dq[newaxis, :]

    # setup unit speed array
    speed = ones([I_dq.shape[0]])

    # reference rotor speed [rpm] and stator currents [A]
    self.input.N0 = speed[0]
    self.input.Id_ref = I_dq[0, 0]
    self.input.Iq_ref = I_dq[0, 1]

    # get the multisimulation
    var_simu = self.get_VarLoadCurrent(I_dq, speed)
    self.var_simu = var_simu

    # add some datakeepers
    Psi_keeper = DataKeeper(
        name="Fluxlinkage",
        symbol="Psi",
        unit="Vs",
        keeper="lambda out: out.mag.Phi_wind_stator.values",
        error_keeper="lambda simu: np.nan",
    )
    var_simu.datakeeper_list.append(Psi_keeper)

    Sd_keeper = DataKeeper(
        name="d-Current Density",
        symbol="Sd",
        unit="A/m**2",
        keeper=f"lambda output: output.elec.Id_ref / {S_wind_act}",
        error_keeper="lambda simu: np.nan",
    )
    var_simu.datakeeper_list.append(Sd_keeper)

    Sq_keeper = DataKeeper(
        name="q-Current Density",
        symbol="Sq",
        unit="A/m**2",
        keeper=f"lambda output: output.elec.Iq_ref / {S_wind_act}",
        error_keeper="lambda simu: np.nan",
    )
    var_simu.datakeeper_list.append(Sq_keeper)

    # run the simulation
    self.run()

    # compute dq flux and store as datakeeper
    output._set_fluxlinkage()

    # restore the machine
    self.machine = org_machine

    return output
