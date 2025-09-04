# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

from abc import ABC, abstractmethod
import numpy as np
import quadprog as mcgibbon_quadprog

class QuadprogInterface(ABC):

    @abstractmethod
    def quadprog(self, H, c, A, b, Aeq, beq):
        """
        abstract method to solve a qp problem
        min 1/2 x^T.H.x + c^T.x
        subj to A.x<=b, Aeq.x=beq

        :param H:
        :param c:
        :param A:
        :param b:
        :param Aeq:
        :param beq:
        :return:
        """
        pass


class McGibbonQuadprog(QuadprogInterface):
    def quadprog(self, H, c, A, b, Aeq, beq):
        # solve_qp  docs:
        # Minimize     1/2 x^T G x - a^T x
        # Subject to   C.T x >= b
        # Parameters
        # ----------
        # G : array, shape=(n, n)
        #    matrix appearing in the quadratic function to be minimized
        # a : array, shape=(n,)
        #    vector appearing in the quadratic function to be minimized
        # C : array, shape=(n, m)
        #    matrix defining the constraints under which we want to minimize the
        #    quadratic function
        # b : array, shape=(m), default=None
        #    vector defining the constraints
        # meq : int, default=0
        #    the first meq constraints are treated as equality constraints,
        #    all further as inequality constraints (defaults to 0).
        # factorized : bool, default=False
        #    If True, then we are passing :math:`R^{−1}` (where :math:`G = R^T R`)
        #    instead of the matrix G in the argument G.

        # Returns
        # -------
        # x : array, shape=(n,)
        #    vector containing the solution of the quadratic programming problem.
        # f : float
        #    the value of the quadratic function at the solution.
        # xu : array, shape=(n,)
        #    vector containing the unconstrained minimizer of the quadratic function
        # iterations : tuple
        #    2-tuple. the first component contains the number of iterations the
        #    algorithm needed, the second indicates how often constraints became
        #    inactive after becoming active first.
        # lagrangian : array, shape=(m,)
        #    vector with the Lagragian at the solution.
        # iact : array
        #    vector with the indices of the active constraints at the solution.

        G_qp = H
        a_qp = -c
        if (Aeq is None) and (A is None):
            x_qp = np.linalg.solve(G_qp, a_qp)
            iact_qp = np.array([])
        else:
            if Aeq is None:
                C_qp = -np.transpose(A)
                b_qp = -b
                meq_qp = 0
            elif A is None:
                C_qp = -np.transpose(Aeq)
                b_qp = -beq
                meq_qp = Aeq.shape[0]
            else:
                C_qp = -np.transpose(np.vstack((Aeq, A)))
                b_qp = - np.hstack((beq, b))
                meq_qp = Aeq.shape[0]
            (x_qp, f_qp, xu_qp, iter_qp, lag_qp, iact_qp) = mcgibbon_quadprog.solve_qp(G_qp, a_qp, C_qp, b_qp, meq_qp)

        return x_qp, iact_qp