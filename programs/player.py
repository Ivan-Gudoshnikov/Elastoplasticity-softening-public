# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import dill
import numpy as np
import matplotlib.pyplot as plt
from solver.springs_view import SpringsView, ColorScheme
from solver.utils import vector_to_matrix
from solver.virial_stress import plot_virial_stress

#A TOOL TO LOAD SAVED SIMULATION DATA FROM .dill, PLAY AND SAVE ANIMATIONS

filename = "large_lattices_examples\\triangles_softening_defect"
#filename_save = "large_lattices_examples\\triangles_softening_defect_stress"

with open(filename+'.dill', 'rb') as handle:
    saved_data = dill.load(handle)

solution = saved_data.get("solution")
process = saved_data.get("process")
viewport = saved_data.get("viewport")

_ = SpringsView(solution, process, viewport, time_text_coords=(-1.1, -0.8))
#_ = SpringsView(solution, process,viewport, time_text_coords=(-1.1, -0.8), filename = filename_save+".mp4", fps=20, color_scheme= ColorScheme.STRESS_ABSOLUTE)#, color_scheme_max_cutoff=600)


## computing the area, assuming rough rectangle:
xi0_mat = vector_to_matrix(solution.XI[:,0], 2)
area = (np.max(xi0_mat[:,0])-np.min(xi0_mat[:,0]))*(np.max(xi0_mat[:,1])-np.min(xi0_mat[:,1]))

figSigma, axSigma = plt.subplots()
axSigma.plot(solution.T, solution.Sigma.T)
axSigma.set(title="Sigma")

figIter, axIter = plt.subplots()
axIter.plot(solution.T[1:], solution.Iters)
axIter.set(title="Number of iterations")
plt.xlabel('t')
plt.ylabel('# of iterations')

plot_virial_stress(process, solution, area)
plt.show()
