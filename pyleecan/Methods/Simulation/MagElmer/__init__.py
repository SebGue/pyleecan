from ....Functions.labels import (
    SHAFTSR_LAB,
    SHAFTSL_LAB,
    YSR_LAB,
    YSL_LAB,
    ROTOR_LAB,
    STATOR_LAB,
    SBS_R_LAB,
    SBS_L_LAB,
    SBR_LAB,
    AS_R_LAB,
    AS_L_LAB,
    AR_LAB,
    AIRBOX_SL_LAB,
    AIRBOX_SR_LAB,
    AIRBOX_R_LAB,
    LAM_LAB,
    YOKE_LAB,
)

R_LAB = ROTOR_LAB + "-0_"
S_LAB = STATOR_LAB + "-0_"
# dictionary to match MagElmer FEA boundary conditions (dict values)
# with pyleecan line boundary properties (dict keys)
# that are set in the build_geometry methods
MagElmer_BP_dict = dict()
MagElmer_BP_dict[R_LAB + AS_R_LAB] = "MASTER_ROTOR_BOUNDARY"
MagElmer_BP_dict[R_LAB + AS_L_LAB] = "SLAVE_ROTOR_BOUNDARY"

MagElmer_BP_dict[R_LAB + SBS_R_LAB] = "MASTER_ROTOR_BOUNDARY"
MagElmer_BP_dict[R_LAB + SBS_L_LAB] = "SLAVE_ROTOR_BOUNDARY"
MagElmer_BP_dict[R_LAB + SBR_LAB] = "SB_ROTOR_BOUNDARY"

MagElmer_BP_dict[R_LAB + YSR_LAB] = "MASTER_ROTOR_BOUNDARY"  # Rotor Yoke Side Right
MagElmer_BP_dict[R_LAB + YSL_LAB] = "SLAVE_ROTOR_BOUNDARY"  # Rotor Yoke Side Left
MagElmer_BP_dict[S_LAB + YSR_LAB] = "MASTER_STATOR_BOUNDARY"  # Stator Yoke Side Right
MagElmer_BP_dict[S_LAB + YSL_LAB] = "SLAVE_STATOR_BOUNDARY"  # Stator Yoke Side Left
MagElmer_BP_dict[SHAFTSR_LAB] = "MASTER_ROTOR_BOUNDARY"  # Shaft Side Right
MagElmer_BP_dict[SHAFTSL_LAB] = "SLAVE_ROTOR_BOUNDARY"  # Shaft Side Left

MagElmer_BP_dict[S_LAB + AS_R_LAB] = "MASTER_STATOR_BOUNDARY"
MagElmer_BP_dict[S_LAB + AS_L_LAB] = "SLAVE_STATOR_BOUNDARY"
MagElmer_BP_dict[S_LAB + AR_LAB] = "AIRGAP_ARC_BOUNDARY"

MagElmer_BP_dict[S_LAB + SBS_R_LAB] = "MASTER_STATOR_BOUNDARY"
MagElmer_BP_dict[S_LAB + SBS_L_LAB] = "SLAVE_STATOR_BOUNDARY"
MagElmer_BP_dict[S_LAB + SBR_LAB] = "SB_STATOR_BOUNDARY"

MagElmer_BP_dict[S_LAB + AIRBOX_SR_LAB] = "MASTER_STATOR_BOUNDARY"
MagElmer_BP_dict[S_LAB + AIRBOX_SL_LAB] = "SLAVE_STATOR_BOUNDARY"

MagElmer_BP_dict[R_LAB + AIRBOX_SR_LAB] = "MASTER_ROTOR_BOUNDARY"
MagElmer_BP_dict[R_LAB + AIRBOX_SL_LAB] = "SLAVE_ROTOR_BOUNDARY"

MagElmer_BP_dict[R_LAB + LAM_LAB + YOKE_LAB] = "ROTOR_EXT_BOUNDARY"
MagElmer_BP_dict[S_LAB + LAM_LAB + YOKE_LAB] = "STATOR_EXT_BOUNDARY"

MagElmer_BP_dict[S_LAB + AIRBOX_R_LAB] = "STATOR_EXT_BOUNDARY"
MagElmer_BP_dict[R_LAB + AIRBOX_R_LAB] = "ROTOR_EXT_BOUNDARY"
