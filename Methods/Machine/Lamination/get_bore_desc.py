from pyleecan.Functions.Machine.Lamination._get_bore_desc import _get_bore_desc

def get_bore_desc(self, sym=1):
    """This method returns an ordered description of the elements 
    that defines the bore radius of the lamination

    Parameters
    ----------
    self : LamSlotMulti
        A LamSlotMulti object

    Returns
    -------
    bore_desc : list 
        list of dictionary with key: "begin_angle", "end_angle", "obj"
    """
    bore_desc = _get_bore_desc(self, sym=sym)

    return bore_desc
