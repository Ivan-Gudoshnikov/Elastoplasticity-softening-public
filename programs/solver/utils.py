# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import numpy as np

def zigzag(period, amplitude, offset):
    def z(t):
        phi=t % period
        if phi<period/4:
            return phi*(4/period)*amplitude + offset
        elif phi<3*period/4:
            return (2-phi*(4/period))*amplitude + offset
        else:
            return (phi*(4/period)-4)*amplitude + offset

    return z

def matrix_to_vector(matrix):
    """
    :param matrix: n x d matrix
    :return: vector ( [d x (1)]  [d x 2] ... [d x (n)])
    """
    (n, d) = matrix.shape
    result = np.zeros(n*d)
    for j in range(0,n):
        result[range(d*j, d*(j+1))] = matrix[j, :]
    return result

def vector_to_matrix(vector, d):
    """
    :param vector: vector ( [d x (1)]  [d x 2] ... [d x (n)])
    :param d:
    :return: matrix: n x d matrix
    """
    n = vector.shape[0] // d
    result = np.zeros((n, d))
    for j in range(0,n):
        result[j, :] = vector[range(d*j, d*(j+1))]
    return result

def tensor_to_matrix(t):
    """
    :param t: m x n x d numpy array
    :return: m x (nd) numpy array: ( t[m x (1) x d]  t[m x (2) x d] ... t[m x (n) x d])
    """
    (m, n, d) = t.shape
    result = np.zeros((m, n * d))
    for j in range(0, n):
        result[:, range(d * j, d * (j + 1))] = t[:, j, :]
    return result


def matrix_to_tensor(matrix, d):
    """
    :param matrix: m x (nd) numpy array: m x (nd) numpy array: ( t[m x (1) x d]  t[m x (2) x d] ... t[m x (n) x d])
    :param d: 3nd dimension size
    :return: m x n x d numpy array
    """
    m = matrix.shape[0]
    n = matrix.shape[1] // d
    result = np.zeros((m, n, d))
    for j in range(0, n):
        result[:, j, :] = matrix[:, range(d * j, d * (j + 1))]
    return result