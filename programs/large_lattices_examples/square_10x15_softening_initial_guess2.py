# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import math
import dill

from solver.elastoplastic_process_linearized import ElastoplasticProcessLinearized
from solver.grid import Grid
import numpy as np
import matplotlib.pyplot as plt

from solver.springs_view import SpringsView

#numbers of nodes along the axes
n1 = 10
n2 = 15
#DEFINITION: logical coordiantes of a node is a tuple (i,j), where i in range(n1), j in range(n2)

def k_func(orig, termin):
    """
    Defines the stiffness (the elasticity parameter) of a spring between two nodes, if such a spring exists
    :param orig: logical coordinates (i,j) of the origin node
    :param termin: logical coordinates (i,j) of the terminus node
    :return: stiffness of the spring
    """
    base_stiffess = 1.
    return base_stiffess

def c0_func(orig, termin):
    """
    Defines the initial yield limit (a plasticity parameter) of a spring between two nodes, if such a spring exists
    :param orig: logical coordinates (i,j) of the origin node
    :param termin: logical coordinates (i,j) of the terminus node
    :return: initial yield of the spring
    """
    base_yeld_stress=0.01
    return base_yeld_stress

def s_func(orig, termin):
    """
    Defines the geometric slope (a plasticity parameter) of a spring between two nodes, if such a spring exists
    :param orig: logical coordinates (i,j) of the origin node
    :param termin: logical coordinates (i,j) of the terminus node
    :return: the geometric slope for the state-dependent flow rule of the spring
    """
    base_s=0.3
    return base_s

def h_func(orig, termin):
    """
    Defines the coefficient of state-dependent feedback (a plasticity parameter) of a spring between two nodes, if such a spring exists
    :param orig: logical coordinates (i,j) of the origin node
    :param termin: logical coordinates (i,j) of the terminus node
    :return: the coefficient of state-dependent feedback for the state-dependent flow rule of the spring
    """
    base_h=1.1
    return base_h

def is_node_func(coords):
    """
    Defines if a particular node should be included in the lattice, allows to make gaps in the n1xn2 lattice
    :param coords: logical coordinates (i,j) of the node
    :return: a boolean value of whether this node should be included
    """
    return True

def xi0_func(coords):
    """
    Defines physical coordinates of a node
    :param coords: logical coordinates (i,j) of the node
    :return: physical coordinates of the node, a touple of real numbers (float)
    """
    x0 = - n1/2.+0.5
    delta = 1.
    (i,j) = coords
    result = (i * delta +x0, j * delta)
    return result

def add_springs_func(orig):
    """
    Defines all of the springs with the origin at a particular node. Allows to specify a spring between particular nodes to exist or not .
    :param orig: logical coordinates (i,j) of the origin node
    :return: a list of logical coordinates of the termina of the springs with such origin
    """
    (i,j) = orig

    termins = []
    if (i < n1-1) and is_node_func((i+1,j)):
        termins.append((i+1,j))
    if (j < n2-1) and is_node_func((i,j+1)):
        termins.append((i,j+1))
    if (i < n1-1) and (j < n2-1) and is_node_func((i+1, j + 1)):
        termins.append((i+1, j + 1))
    if (i > 0) and (j < n2 - 1) and is_node_func((i - 1, j + 1)):
        termins.append((i - 1, j + 1))
    return termins


def r_and_r_prime_component_functions(node_coords, component):
    """
    Defines the displacement boundary condition at a node
    :param node_coords: logical coordinates (i,j) of the node
    :param component: 0 is x, 1 is y; allows to turn on/off constraints on individual components
    :return: Function of time, representing the component of the displacement and its derivative at the node
    """
    func = None
    if node_coords[1] == 0: #lower row of nodes
        if component == 0: #for x component
            func = None #free motion along x axis
        else: #for y component
            func = lambda t: (0,0) # hold
    elif node_coords[1] == n2-1: #upper row of nodes
        if component == 0: #for x component
            func = None #free motion along x axis
        else: #for y component
            rate = 0.1 #pull up
            func = lambda t: (rate*t, rate)

    if node_coords[0] == math.ceil(n1/2.)-1: #for the nodes on the vertical axis of symmetry
        if (node_coords[1] == 0): # lower row
            if component == 0: #for x component
                func = lambda t: (0,0) # hold

    return func

def force_func(node_coords):
    """
    Defines an external force at a node
    :param node_coords: logical coordinates (i,j) of the node
    :return: Vector function (touple of floats), representing the components of the external force at the node
    """
    force = lambda t: (0,0)
    return force

grid = Grid(n1, n2, is_node_func, xi0_func, add_springs_func, k_func, c0_func, r_and_r_prime_component_functions, force_func, s_func, h_func)
process = ElastoplasticProcessLinearized(grid.Q, grid.xi, grid.k, grid.c0, grid.s, grid.h, grid.d,
                                         grid.boundary_condition.R, grid.boundary_condition.r, grid.boundary_condition.f, True)



t0 = 0
dt = 0.01
nsteps = 1000
e0=grid.e0
a0=grid.a0
p0=grid.p0
xi0=grid.xi

#using initial guess with the special offset
v = np.zeros(process.m)
for w in range(n1):
    node_id_1  = grid.node_id_by_coord[w, n2-4-w]
    node_id_2 = grid.node_id_by_coord[w, n2 - 3 - w]
    spring_id = grid.connections.index((node_id_1,node_id_2))
    v[spring_id] = 2e-5

def initial_guess_a_offset(i):
    if i>275 and i<285:
        return v
    else:
        return np.zeros(process.m)

solution = process.solve_iteratively_in_V(e0, t0, dt, nsteps, a0, p0)

figSigma, axSigma = plt.subplots()
axSigma.plot(solution.T, solution.Sigma.T)
axSigma.set(title="Sigma")

figA, axA = plt.subplots()
axA.plot(solution.T, solution.A.T)
axA.set(title="A")

figIter, axIter = plt.subplots()
axIter.plot(solution.T[1:], solution.Iters)
axIter.set(title="Number of iterations")
plt.xlabel('t')
plt.ylabel('# of iterations')


viewport = ((-5.5, 5.5), (-1, 15))
filename = "square_10x15_softening_initial_guess2"

with open(filename+'.dill', 'wb') as handle:
    dill.dump({"grid": grid, "process": process, "solution": solution, "viewport": viewport}, handle, protocol=dill.HIGHEST_PROTOCOL)

sv = SpringsView(solution, process, viewport, time_text_coords=(-1.1, -0.8))
#USE THE NEXT ONE INSTEAD TO SAVE THE ANIMATION INTO .mp4, REQUIRES FFMPEG INSTALLED IN THE SYSTEM
#sv = SpringsView(solution, process,viewport, time_text_coords=(-1.1, -0.8), filename=filename+".mp4",fps=20)

plt.show()
