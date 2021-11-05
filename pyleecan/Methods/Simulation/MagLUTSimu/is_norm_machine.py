def is_norm_machine(self, machine):
    """Method to check if a machine is equal to the normlized simulation machine."""
    norm_input_machine = self.get_norm_machine(machine)
    norm_input_machine.stator.winding = None

    # get the norm. machine
    simu_machine = self.get_norm_machine(self.machine)
    simu_machine.stator.winding = None

    diff_list = simu_machine.compare(norm_input_machine)
    filter = ["logger_name", "frame", "shaft", "path", "name", "desc"]
    diff_filtered = []
    for diff in diff_list:
        if diff.split(".")[-1].split(" ")[0] not in filter:
            diff_filtered.append(diff)

    if not diff_filtered:
        return True

    return False
