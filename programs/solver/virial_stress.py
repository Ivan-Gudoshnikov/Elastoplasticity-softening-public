# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import numpy as np
import matplotlib.pyplot as plt
from solver.elastoplastic_process_linearized import matrix_of_directions_D, phi, ElastoplasticProcessLinearized

def get_virial_stress(process, area, Sigma):
    D = matrix_of_directions_D(process.Q, process.xi0, process.d)
    phi0 = phi(process.Q, process.xi0, process.d)

    virial_stress = np.zeros((process.d, process.d, Sigma.shape[1]))

    for i in range(Sigma.shape[1]):
        vs = D.T @ np.diag(Sigma[:, i]) @ np.diag(phi0) @ D / area
        virial_stress[:, :, i] = vs[:, :]

    return virial_stress


def plot_virial_stress(process: ElastoplasticProcessLinearized, solution: ElastoplasticProcessLinearized.Solution, area, ax = None):
    """
    Plot the 2,2-component of the total stress, formula (112)

    :param process:
    :param solution:
    :param area:
    :param ax:
    :return:
    """
    if ax is None:
        _, axSigma_total = plt.subplots()
        axSigma_total.set(title="Total stress along the vertical direction")
        plt.xlabel('t')
        plt.ylabel('\sigma_{total}^{22}')
    else:
        axSigma_total = ax
        plt.sca(ax)

    virial_stress_mat = get_virial_stress(process, area, solution.Sigma)
    virial_stress_vec = virial_stress_mat.reshape(-1, virial_stress_mat.shape[-1])

    #axSigma_total.plot(solution.T, virial_stress_vec[[0,1,3],:].T, label=[r'$\sigma^{11}_{total}$',r'$\sigma^{12}_{total}=\sigma^{21}_{total}$', r'$\sigma^{22}_{total}$'])
    axSigma_total.plot(solution.T, virial_stress_vec[3, :].T)
    #axSigma_total.legend(loc='center right')

    return axSigma_total




