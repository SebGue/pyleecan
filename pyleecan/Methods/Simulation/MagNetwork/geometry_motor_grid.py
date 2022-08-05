# -*- coding: utf-8 -*-

# from tkinter.ttk import _TreeviewColumnDict
from tkinter import *
from turtle import st
import numpy as np
from pyleecan.Functions.load import load
from scipy.fft import ifftn
import matplotlib.pyplot as plt
from collections import defaultdict
import re
import math

#######################################################################################
# Store the motor points coords (from ".geo") in a dict and convert them to radial
#######################################################################################
def lines_that_start_with(string, input_file):
    """returns the list of elements of a line (from a file) starting by a string
    Parameters
    ----------
    string : str
        The string to be searched for
    input_file : str
        The file path of the ".geo"

    Returns
    -------
    list_point_file : list
        list of coords defined in point(i)
    """
    for line in input_file:
        if line.startswith(string):
            # return re.split(r"[{,\s]", line.strip())
            list_point_file = re.split(r"[{,\s]", line.strip())
            return list_point_file


def count_word_in_file(file, word):
    """counts the number of a word in a file
    Parameters
    ----------
    file : str
        The file path
    word : str
        The string to be searched for

    Returns
    -------
    list_point_file : list
        list of coords defined in point(i)
    """
    datafile_data = file.read()
    word_count = datafile_data.count(word)
    return word_count


def geometry_points_cartesian(geometry_file_path, count_word_occurence):
    """returns the position of a point defined by an r and a theta in the grid
    Parameters
    ----------
    geometry_file_path : str
        The file path of the ".geo"
    count_word_occurence : integer
        The number of occurence of a string

    Returns
    -------
    coords : dict
        dict of all points coords (x,y)
    """
    file = open(geometry_file_path, "r")
    coords = {}
    for i in range(count_word_occurence):
        list_point_geometry = lines_that_start_with("Point(" + str(i + 1) + ")", file)
        list_point_geometry_x = float(list_point_geometry[3])
        list_point_geometry_y = float(list_point_geometry[5])
        coords["pt" + str(i + 1)] = [list_point_geometry_x, list_point_geometry_y]
    return coords


def convert_cartesian_to_radial(x, y):
    """converts the cartesian coords of a point (x,y) to radial(r,theta)
    Parameters
    ----------
    x : float
        The x-ccordinate of the point
    y : float
        The y-ccordinate of the point

    Returns
    -------
    theta : float
        The theta-ccordinate of the point
    r : float
        The r-ccordinate of the point
    """
    theta = math.atan(y / x) * (180 / np.pi)  # Theta in deg
    r = (x * 1000) / math.cos(theta / (180 / np.pi))  # r in mm
    return theta, r


def geometry_points_polar(cartesian_coords):
    """converts the cartesian coords of a point (x,y) to radial(r,theta)
    Parameters
    ----------
    cartesian_coords : dict
        The dictionary containing points with cartesian coords (x,y)

    Returns
    -------
    polar_coords : dict
        The dictionary containing points with cartesian coords (x,y)
    """
    polar_coords = {}
    for ll in range(len(cartesian_coords)):
        polar_coords["pt" + str(ll + 1)] = list(
            convert_cartesian_to_radial(
                cartesian_coords["pt" + str(ll + 1)][0],
                cartesian_coords["pt" + str(ll + 1)][1],
            )
        )
    return polar_coords


#######################################################################################
# Definition of the position of a point from the motor geometry in the grid
#######################################################################################
def position_point(machine_coords, dtheta, dr, R_shaft):
    """returns the position of a point defined by an r and a theta in the grid
    Parameters
    ----------
    machine_coords : dict
        Dictionary containing the points poisition in the geometry/ point = [theta_machine, r_machine]
    dtheta : integer
        Disctretization step according to the theta-axis
    dr : integer
        Discretization step according to the r-axis

    Returns
    -------
    grid_coords : dict
        Dictionary containing the points position in the grid/ point = [theta_grid, r_grid]
    """

    # Initialization of position_grid
    grid_coords = {}

    # Computing theta_grid and r_grid from theta_machine and r_machine
    for oo in range(len(machine_coords)):

        # Theta from machine geometry to grid
        theta_machine = machine_coords["pt" + str(oo + 1)][0]
        theta_grid = (round(theta_machine / dtheta) * dtheta) / dtheta

        # r from machine geometry to grid
        r_machine = machine_coords["pt" + str(oo + 1)][1]
        r_grid = (round(r_machine / dr) * dr - R_shaft) / dr

        # Storing points in a dict
        grid_coords["pt" + str(oo + 1)] = [theta_grid, r_grid]

    return grid_coords


#######################################################################################
# Definition of the electric motor geometry : "geometry_motor method"
#######################################################################################
def geometry_motor(self, geometry_file_path, dtheta, dr):
    """Defines the discretized geometry of the electric machine under study

    Parameters
    ----------
    ..

    Returns
    -------
    ..
    """
    #######################################################################################
    # Definition of the general machine input parameters
    #######################################################################################

    # Getting the machine object
    Machine = self.parent.machine

    # Conversions
    rad_to_deg = 180 / np.pi
    m_to_mm = 1000

    # Machine periodicity
    if Machine.comp_periodicity_spatial()[1] == True:
        periodicity = Machine.comp_periodicity_spatial()[0]
    else:
        periodicity = Machine.comp_periodicity_spatial()[0] / 2

    # Pole pitch angle
    angle_tp = (np.pi / periodicity) * rad_to_deg

    # Number of PMs per period
    nb_PM_per_period = round(Machine.rotor.get_pole_pair_number() / periodicity)

    # Number of stator teeth per periods
    nb_stator_teeth_per_period = round(Machine.stator.get_Zs() / (2 * periodicity))

    # Number of winding layers
    nb_layers = Machine.stator.winding.Nlayer

    # Active length of the motor, assuming that la(rotor) = la(stator)
    la = Machine.rotor.L1

    #######################################################################################
    # Definition of the angles of the motor elements in "deg"
    #######################################################################################

    """ Angles defining the stator elements """
    # Angular opening of the stator slot
    angle_slot = nb_layers * (
        Machine.stator.slot.comp_angle_active_eq() * rad_to_deg / nb_layers
    )

    # Angle of a half stator tooth
    angle_half_tooth = (angle_tp - nb_stator_teeth_per_period * angle_slot) / (
        2 * nb_stator_teeth_per_period
    )

    # Angle of a stator tooth
    angle_tooth = 2 * angle_half_tooth

    """ Angles defining the rotor elements """
    # Angular opening of the rotor magnet
    angle_magnet = Machine.rotor.slot.comp_angle_active_eq() * rad_to_deg

    # Angle of half the air-gap-PM : angle of one half airgap between 2 adjascent PMs
    angle_half_airgap_PM = (angle_tp - angle_magnet * nb_PM_per_period) / (
        2 * nb_PM_per_period
    )

    # Angle of the air-gap-PM : angle of the airgap between 2 adjascent PMs
    angle_airgap_PM = 2 * angle_half_airgap_PM

    #######################################################################################
    # Definition of the radius of the motor elements in "mm"
    #######################################################################################

    # Rotor interior and exterior radius
    radius_rotor_interior = Machine.rotor.Rint * m_to_mm
    radius_rotor_exterior = (
        radius_rotor_interior + Machine.rotor.comp_height_yoke() * m_to_mm
    )  # rotor outer radius without magnets
    radius_PM = (
        Machine.rotor.comp_radius_mec() * m_to_mm
    )  # exterior radius including the PMs

    # Stator interior and exterior radius
    radius_stator_interior = Machine.stator.Rint * m_to_mm
    radius_stator_exterior = Machine.stator.Rext * m_to_mm

    # Stator slots interior and exterior radius
    radius_stator_slot_exterior = (
        radius_stator_exterior - Machine.stator.comp_height_yoke() * m_to_mm
    )

    radius_stator_slot_interior = (
        radius_stator_slot_exterior - Machine.stator.slot.comp_height_active() * m_to_mm
    )

    #######################################################################################
    # Defining the geometry points in the position_machine dictionary
    #######################################################################################

    # Opening the file from the geometry file path
    file = open(geometry_file_path, "r")

    # Defining the number of ocurrence of "Point" in the file path
    count_point_occurence = count_word_in_file(file, "Point")

    # Getting the Cartesian points coords in a dict "cartesian_coords"
    cartesian_coords = geometry_points_cartesian(
        geometry_file_path, count_point_occurence
    )

    # Defining the radial points coords from the cartesian ones in a dict "geometry_polar_coord"
    geometry_polar_coord = geometry_points_polar(cartesian_coords)

    # Defining the position of points in the grid
    grid_polar_coords = position_point(
        geometry_polar_coord, dtheta, dr, radius_rotor_interior
    )

    #######################################################################################
    # Plotting the position of the points in the geometry for verification
    #######################################################################################
    fig = plt.figure()
    ax = fig.add_subplot(projection="polar")

    for jj in range(len(grid_polar_coords)):
        ax.scatter(
            grid_polar_coords["pt" + str(jj + 1)][0] * (np.pi / 180),  # Theta
            grid_polar_coords["pt" + str(jj + 1)][1],  # r
            s=15,  # width of the point
        )
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    plt.xlabel("r position")
    plt.show()

    #######################################################################################
    # Grid Definition
    #######################################################################################

    # Definition of the number of elements theta and r
    N_element_theta = round(angle_tp / dtheta)
    N_element_r = round((Machine.stator.Rext - Machine.rotor.Rint) / dr)
    N_total_element = N_element_theta * N_element_r

    # Definition of the r and theta axes of the grid
    theta_axis = np.linspace(0, (N_element_theta + 1), dtheta)
    r_axis = np.linspace(0, (N_element_r + 1), dr)

    # Plot the grid of size (N_element_theta * N_element_r)
    plt.title("Discretization grid")
    plt.xlabel("N_element_theta")
    plt.ylabel("N_element_r")

    plt.plot(
        theta_axis,
        r_axis,
    )

    plt.clf()  # clear figure
    # Define the x an y axes steps of the grid
    plt.xticks(np.arange(min(theta_axis), max(theta_axis) + 1, dtheta))
    plt.yticks(np.arange(min(r_axis), max(r_axis) + 1, dr))

    plt.grid()
    plt.show()

    # Attribute an array of N_total_element size of zeros to cells_materials
    cells_materials = np.zeros(N_total_element, dtype=np.float64)

    return grid_polar_coords
