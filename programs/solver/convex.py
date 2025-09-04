# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).


from abc import ABC, abstractmethod
import numpy as np
import scipy.optimize
import warnings

from solver.quadprog_interface import QuadprogInterface, McGibbonQuadprog

#warnings.filterwarnings('error')
class ConvexSet(ABC):

    @abstractmethod
    def projection(self, H, x, method):
        """
        abstract method for finding a closest point to a set
        :param H: inner product x^T.H.y  defines the norm
        :param x: the point to project
        :param method: link to a respective method to use
        :return: the closest point
        """
        pass

    @abstractmethod
    def __contains__(self, x):
        """
        Implement a check if a point x is from the set
        :param x:
        :return:
        """
        pass

    @abstractmethod
    def normal_cone(self, H, x):
        """
        the normal cone to the set at x in the sense of the inner product  x^T.H.y
        :param H: inner product x^T.H.y  needed to construct the normal cone
        :param x:
        :return: ConvexCone
        """
        pass

    @abstractmethod
    def tangent_cone(self, x):
        """
        :param x:
        :return:
        """

class Polytope(ConvexSet):
    eps = 1e-12

    def __init__(self, A, b, Aeq, beq):
        """
        Create a polytope in R^n defined via inequality constraints A @ x <= b and equality constraints Aeq @ x = beq
        :param A:
        :param b:
        :param Aeq:
        :param beq:
        """
        if A is not None:
            self.n = A.shape[1]
        elif Aeq is not None:
            self.n = Aeq.shape[1]
        else:
            self.n = None  # Means any :-(
        self.A = A
        self.b = b
        self.Aeq = Aeq
        self.beq = beq

    def projection(self, H, x, method: QuadprogInterface):
        """
        method for finding a closest point to a set using given implementation of a quadratic programming algorithm
        :param H: inner product x^T.H.y  defines the norm
        :param x: the point to project
        :param method: link to a implementation of a quadprog algorithm
        :return: the closest point
        """

        # (proj-x)^T.H.(proj-x) = proj^T.H.proj - 2.x^T.H.proj + x^T.H.x
        # 1/2 proj^T.H.proj + (- H.x)^T.proj
        #  1/2 x^T.H.x + c^T.x

        if (self.A is None) and (self.Aeq is None):
            proj = x
        else:
            proj, active = method.quadprog(H, -H @ x, self.A, self.b, self.Aeq, self.beq)
        return proj

    def linprog(self, c):
        result = scipy.optimize.linprog(c, self.A, self.b, self.Aeq, self.beq)
        if not result.success:
            raise NameError("Can't do linprog, errorcode", result.status)
        return result

    def __contains__(self, x):
        #return self.contains_via_quadprog_projection(x)
        return self.contains_direct(x)

    def contains_direct(self, x):
        if x.shape[0] != self.n:
            raise NameError("x is of wrong dimension!")
        if self.A is not None:
            diff = self.A @ x - self.b
            ineq = all(diff <= self.eps)
        else:
            ineq = True
        if self.Aeq is not None:
            diff_eq = self.Aeq @ x - self.beq
            eq = (np.linalg.norm(diff_eq, ord=np.inf) < self.eps)
        else:
            eq = True
        return ineq and eq

    def contains_via_quadprog_projection(self,x):
        x1 = self.projection(np.eye(self.n), x, McGibbonQuadprog())
        return np.linalg.norm(x-x1)<self.eps



    def normal_cone(self, H, x):
        raise NameError(
            "NOT IMPLEMENTED YET: To find the constraints of the normal cone to a polytope we need to solve a constraint enumeration problem ")

    def tangent_cone(self, x):
        """
        :param x:
        :return:
        """
        if x not in self:
            raise NameError("x is not from the set!")
        cone_Aeq = self.Aeq
        if self.Aeq is not None:
            cone_beq = np.zeros_like(self.beq)
        else:
            cone_beq = None

        if self.A is not None:
            ineq = np.abs(np.matmul(self.A, x) - self.b) <= self.eps
            if np.any(ineq):
                cone_A = self.A[ineq, :]
                cone_b = np.zeros(np.sum(ineq))
            else:
                cone_A = None
                cone_b = None
        else:
            cone_A = None
            cone_b = None
        return Polytope(cone_A, cone_b, cone_Aeq, cone_beq)
