# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from solver.elastoplastic_process_linearized import ElastoplasticProcessLinearized
from solver.quadprog_interface import McGibbonQuadprog

def n_func(lambdas_plus, lambdas_minus, process:ElastoplasticProcessLinearized):
    return process.plasticity_data.normals_minus @ lambdas_minus + process.plasticity_data.normals_plus @ lambdas_plus



def example0_projection_chart(process:ElastoplasticProcessLinearized, solution,i_to_examine, width, height, arrow_base_w=None, draw_next_step=True, draw_iterations=True):
    eps = 1e-10
    if arrow_base_w is None:
        arrow_base_w = width/60.

    t_0 = solution.T[i_to_examine]
    yva_0 = solution.YvA[:, i_to_examine]
    a_0 = yva_0[process.dimV:]
    print("projection from: ", a_0)
    t_1 = solution.T[i_to_examine + 1]

    if draw_next_step:
        yva_1r = solution.YvA[:, i_to_examine + 1]
        a_1r = yva_1r[process.dimV:]
        print("fixed point: ", a_1r)

    k=16
    X = np.linspace(a_0[0]-0.1*width, a_0[0]+1.1*width, k)
    Y = np.linspace(a_0[1]-0.1*height, a_0[1]+1.1*height, k)
    U, V = np.meshgrid(X, Y)
    U1 = np.zeros_like(U)
    V1 = np.zeros_like(V)
    quadprog_solver = McGibbonQuadprog()

    figProjChart, axProjChart = plt.subplots()
    axProjChart.set(title="Projections", aspect='equal', xlim=(X[0], X[-1]),ylim=(Y[0], Y[-1]))

    for i in range(k):
        for j in range(k):
            a_1=np.array([U[i, j], V[i,j]])
            C_1 = process.moving_set_func(t_1, a_1)
            yva_1proj = C_1.projection(process.S_V, yva_0, quadprog_solver)
            a_1proj = yva_1proj[process.dimV:]
            U1[i, j] = a_1proj[0]
            V1[i, j] = a_1proj[1]

            vec = a_1-a_1proj
            if np.linalg.norm(vec)<eps:
                axProjChart.add_line(Line2D([U[i,j], U[i,j]], [V[i,j], V[i,j]], marker='o',
                              markerfacecolor='c', markeredgecolor='c', linestyle='None'))
            else:
                #axProjChart.add_line(Line2D([U[i,j], U1[i,j]], [V[i,j], V1[i,j]]))
                #plt.arrow(U[i, j], V[i,j],U1[i, j]-U[i, j], V1[i, j]-V[i,j],
                #          length_includes_head=True, head_width=1e-6,head_length=2e-6)
                vec_ort = np.array([vec[1], -vec[0]])/np.linalg.norm(vec)
                xy = np.array([a_1proj, a_1+0.5*arrow_base_w*vec_ort,a_1-0.5*arrow_base_w*vec_ort])
                arrow_color = (0.12, 0.47, 0.71, 0.7)

                axProjChart.add_patch(plt.Polygon(xy, fc=arrow_color))

    #axProjChart.quiver(X, Y, U1, V1)


    if draw_next_step:
        solution_markers = Line2D([a_1r[0], a_1r[0]], [a_1r[1], a_1r[1]], marker='o',
                              markerfacecolor='r', markeredgecolor='r', linestyle='None')

    if draw_iterations:
        iters_xdata = solution.iteration_steps[i_to_examine][1, :]
        iters_ydata = solution.iteration_steps[i_to_examine][2, :]
    else:
        #draw the initial point only
        iters_xdata = [a_0[0],a_0[0]]
        iters_ydata = [a_0[1],a_0[1]]

    iters_line = Line2D(iters_xdata, iters_ydata, marker='o',
                        markerfacecolor='k', markeredgecolor='k', linestyle='-', color='k')
    axProjChart.add_line(iters_line)

    if draw_next_step:
        axProjChart.add_line(solution_markers)

    plt.xlabel('a_1')
    plt.ylabel('a_2')

    return axProjChart





