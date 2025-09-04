# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

from enum import Enum
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.animation as animation

from solver.elastoplastic_process_linearized import ElastoplasticProcessLinearized
from solver.utils import vector_to_matrix

class ColorScheme(Enum):
    STRESS_TO_YIELD = 0
    STRESS_ABSOLUTE = 1
    DAMAGE = 2
    TOTAL_ELONGATION = 3
    PLASTIC_ELONGATION = 4

class SpringsView:
    """
    Everything to plot the animated latices
    """

    eps = 1e-12

    def __init__(self, solution, problem:ElastoplasticProcessLinearized, lim, time_text_coords, time_text_format='.2f', filename=None, fps=None, dpi=300, y_offset_for_1d=0, color_scheme:ColorScheme = ColorScheme.DAMAGE, color_scheme_max_cutoff=-1):
        """
        @param solution:solution to draw
        @param problem:
        @param lim:
        @param time_text_coords:
        @param filename:
        @param fps:
        """
        if (problem.d != 1) and (problem.d != 2):
            raise NameError("3d networks are not supported yet")
        self.solution = solution
        self.problem = problem
        self.color_scheme = color_scheme
        self.color_scheme_max_cutoff = color_scheme_max_cutoff
        self.fig, self.ax = plt.subplots()
        # self.ax.grid()
        self.ax.set(xlim=lim[0], ylim=lim[1], aspect='equal')

        self.nodes_markers = Line2D([0, 1], [0, 1], marker="None", linestyle="None", markerfacecolor="k",
                                    markeredgecolor="k", markersize=5)
        self.ax.add_line(self.nodes_markers)

        self.springs_lines = []
        for i in range(problem.m):
            self.springs_lines.append(Line2D([0, 1], [0, 1], marker=None, color="k", linewidth=1))
            self.ax.add_line(self.springs_lines[-1])  # add the last created line to the axes

        self.artists = self.springs_lines.copy()
        self.artists.append(self.nodes_markers)

        self.time_text = plt.text(time_text_coords[0], time_text_coords[1], "t=")
        self.time_text_format=time_text_format
        self.artists.append(self.time_text)

        self.max_Sigma = np.max(np.abs(self.solution.Sigma[:, :self.color_scheme_max_cutoff]))


        self.max_A = np.max(np.abs(self.solution.A[:, :self.color_scheme_max_cutoff]))
        self.max_X = np.max(np.abs(self.solution.X[:, :self.color_scheme_max_cutoff]))
        self.max_P = np.max(np.abs(self.solution.P[:, :self.color_scheme_max_cutoff]))

        def init_animation():
            return self.artists

        def update_animation(i):
            self.time_text.set_text("t=" + format(self.solution.T[i], self.time_text_format))
            if problem.d == 2:
                xi = vector_to_matrix(self.solution.XI[:, i], self.problem.d)
            elif problem.d == 1:
                xi = np.vstack((self.solution.XI[:, i], np.zeros_like(self.solution.XI[:, i]))).T

            self.nodes_markers.set_data(xi[:, 0], xi[:, 1]) #to impmement: y_offset_for_1d

            for j in range(problem.m):
                xdata = [xi[self.problem.connections[j][0], 0], xi[self.problem.connections[j][1], 0]]
                ydata = [xi[self.problem.connections[j][0], 1]+y_offset_for_1d*j, xi[self.problem.connections[j][1], 1]+ y_offset_for_1d*j]
                self.springs_lines[j].set_data(xdata, ydata)

                (hue, thickness) = self.values_to_update_linear_softening(i, j)

                self.springs_lines[j].set_color(matplotlib.colors.hsv_to_rgb((abs(hue), 1, 0.9)))
                self.springs_lines[j].set_linewidth(thickness)
            return self.artists

        self.ani = animation.FuncAnimation(self.fig, update_animation, init_func=init_animation, frames=self.solution.T.shape[0],
                                           interval=1, blit=True, repeat=True)
        if filename is not None:
            self.ani.save(filename, fps=fps,dpi = dpi)

    def values_to_update_linear_softening(self,i,j):
        upper_sigma= self.problem.c0[j] + (1-self.problem.h[j])*self.problem.s[j]*self.solution.A[j, i]
        lower_sigma=-upper_sigma
        active = 0
        if abs(self.solution.Sigma[j, i] - upper_sigma) < self.eps:
            active = 1
        elif abs(self.solution.Sigma[j, i] -lower_sigma) < self.eps:
            active = -1

        if active == 0:
            thickness = 2
        else:
            thickness = 4

        match self.color_scheme:
            case ColorScheme.STRESS_TO_YIELD:
                hue1 = abs(self.solution.Sigma[j, i]) / upper_sigma
            case ColorScheme.STRESS_ABSOLUTE:
                hue1 = abs(self.solution.Sigma[j, i]) / self.max_Sigma
            case ColorScheme.DAMAGE:
                hue1 = abs(self.solution.A[j, i]) / self.max_A
            case ColorScheme.TOTAL_ELONGATION:
                hue1 = abs(self.solution.X[j, i]) / self.max_X
            case ColorScheme.PLASTIC_ELONGATION:
                hue1 = abs(self.solution.P[j, i]) / self.max_P

        hue2 = (1 - min([max([-1, hue1]), 1])) * 0.25  # linear interpolation between green(0.25) and red (0)
        return hue2, thickness


