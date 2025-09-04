# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import sys
import numpy as np
import scipy
from solver.convex import Polytope
from solver.quadprog_interface import McGibbonQuadprog

def phi(Q, xi, d):
    """
    Computes lengths of the springs at configuration xi, see formula (42) in the paper
    :param Q: incidence matrix
    :param xi: positions of the nodes, a (n*d)-vector, xi[d*j + k]$ is the k-th coordinate of the node j, k in range(d), j in range(n)
    :param d: the number of spatial dimensions
    :return: a m-vector of lengths of the springs
    """
    n = Q.shape[0]
    m = Q.shape[1]
    result = np.zeros(m)
    J = range(n)
    K = range(d)
    for i in range(m):
        result[i] = np.sqrt(np.sum([np.square(np.dot(xi[[d * j + k for j in J]], Q[:, i])) for k in K]))
    return result

def matrix_of_directions_D(Q, xi0, d):
    """
    A m x d matrix of unit vectors in the directions from the terminus to the origin of the springs arranged as rows, see formula (44) in the paper
    :param Q: incidence matrix
    :param xi0: positions of the nodes, a (n*d)-vector, xi[d*j + k]$ is the k-th coordinate of the node j, k in range(d), j in range(n)
    :param d: the number of spatial dimensions
    :return: the m x d matrix
    """
    phi0 = phi(Q,xi0,d)
    n = Q.shape[0]
    m = Q.shape[1]
    result = np.zeros((m,d))
    J = range(n)
    for i in range(m):
        for k in range(d):
            result[i, k] = np.dot(xi0[[d * j1 + k for j1 in J]], Q[:, i])/phi0[i]
    return result

def d_xi_phi(Q, xi0, d):
    """
    Jacobi matrix of the function xi at a reference configuration xi0, see formula (43) in the paper
    :param Q: incidence matrix
    :param xi0: reference configuration to linearize at - positions of the nodes, a (n*d)-vector, xi0[d*j + k]$ is the k-th coordinate of the node j, k in range(d), j in range(n)
    :param d: the number of spatial dimensions
    :return: m x nd Jacobi matrix
    """
    n = Q.shape[0]
    m = Q.shape[1]
    phi0 = phi(Q, xi0, d)
    result = np.zeros((m, n * d))
    J = range(n)
    K = range(d)
    for i in range(m):
        for j in J:
            for k in K:
                result[i, d * j + k] = np.dot(xi0[[d * j1 + k for j1 in J]], Q[:, i]) * Q[j, i] / phi0[i]
    return result

class ElastoplasticProcessLinearized:
    class Solution:
        pass

    def __init__(self, Q, xi0, k, c0, s, h, d, R, r, f, demand_enough_constraints=True):
        """
        :param Q: incidence matrix
        :param xi0: reference configuraton
        :param k: vector of stiffnesses (elasticity parameters)
        :param c0: vector of initial yield limits (plasticity parameters)
        :param s: vector of the geometric slope (plasticity parameters)
        :param h: vector of the coefficients of state-dependent feedback (plasticity parameters)
        :param d: amount of spatial dimensions
        :param R: matrix from additional constraint R(xi0+zeta) + r(t)=0
        :param r: function r(t) from the additional constraint
        :param f: external force applied at the nodes - function of t
        :param demand_enough_constraints: check the kinematic determinacy
        """

        self.Q = Q

        self.k = k
        self.K = np.diag(k)
        self.Kinv = np.diag(1. / k)

        self.c0 = c0
        self.s = s
        self.h = h
        self.S = np.diag(s)
        self.Sinv = np.diag(1. / s)
        self.H = np.diag(h)

        # amount of springs
        self.m = Q.shape[1]
        # amount of nodes
        self.n = Q.shape[0]
        # spatial dimension
        self.d = d
        # amount of additional constraints
        self.q = R.shape[0]

        if R.shape[1] != self.n * self.d:
            raise NameError("Wrong matrix dimension of R !")
        if np.linalg.matrix_rank(R) != self.q:
            raise NameError("Dependent constraints in R !")

        self.R = R
        self.r = r
        self.f = f

        # explicit list of edges with endpoints, computed from Q
        self.connections = []
        for i in range(self.m):
            spring_tuple = (np.where(self.Q[:, i] == 1)[0][0], np.where(self.Q[:, i] == -1)[0][0])
            self.connections.append(spring_tuple)

        self.xi0 = xi0

        self.d_xi_phi = d_xi_phi(Q, xi0, d)

        self.R0 = scipy.linalg.null_space(self.R)
        self.Rp = np.linalg.pinv(self.R)

        # computing basis in \mathcal{U}
        ubasis_candidate = self.K @ self.d_xi_phi @ self.R0 #see formula (49)
        if demand_enough_constraints:
            constraints_deficit = ubasis_candidate.shape[1] - np.linalg.matrix_rank(ubasis_candidate)
            if constraints_deficit > 0:  # if columns of result are not linearly independent, throw an error
                raise NameError("Not enough constraints for kinematic determinacy! Deficit:" + str(constraints_deficit))
        # ubasis = scipy.linalg.orth(ubasis_candidate) #orthogonalize, if necessary
        self.ubasis = ubasis_candidate

        # To compute F1 take the "upper" part of the pseudo-inverse matrix of the combined  matrix [(d_xi_phi)^T  R^T], see formula (48)
        self.full_EM = np.hstack((self.d_xi_phi.T, self.R.T))  # equilibrium matrix of the constrained system
        self.full_EM_pinv = np.linalg.pinv(self.full_EM)
        self.F1 = self.full_EM_pinv[0:self.m, :]

        self.full_CM_pinv = self.full_EM_pinv.T

        #computing basis in \mathcal{V}, see Proposition 4.1 i)
        self.vbasis = scipy.linalg.null_space((self.Kinv @ self.ubasis).T)
        #dimensions of the fundamental subspaces
        self.dimU = self.ubasis.shape[1]
        self.dimV = self.vbasis.shape[1]

        if (self.dimU + self.dimV != self.m) and demand_enough_constraints:
            raise NameError("Dimensions of U and V do not add up!")

        # COMPUTING P_U, P_V in coodinates of U and V:
        #a way to compute the projection in terms of coordinates is to invert the transformation and ditch the components from the other subspace:
        stacked_inv = np.linalg.inv(np.hstack((self.ubasis, self.vbasis)))
        self.P_U_coords = stacked_inv[range(0, self.dimU), :]
        self.P_V_coords = stacked_inv[range(self.dimU, self.m), :]

        # pre-computing values for the sweeping process
        self.G_V = self.P_V_coords @ self.K @ self.d_xi_phi @ self.Rp
        self.G = self.vbasis @ self.G_V

        self.F = self.ubasis @ self.P_U_coords @ self.F1


        print("Pre-computations are complete.")
        print("m = ", self.m, ", n = ", self.n, ", d = ", self.d, ", q = ", self.q)
        print("dim U = ", self.dimU, " dim V = ", self.dimV)

    def solve_iteratively_in_V(self, e0, t0, dt, nsteps, a0, p0, save_iterations=False, initial_guess_a_offset=None,  iter_max = 5000, value_max = 1e+12, eps = 1e-12):
        """
        Implicit catch-up scheme of reduced dimension, Algorithm 2
        :param e0: initial values for elastic elongation, m-vector
        :param t0: initial time
        :param dt: time-step
        :param nsteps: number of time-steps
        :param a0: initial values for the damage variable, m-vector
        :param p0: initial values for plastic elongation, m-vector
        :param save_iterations: set to True to include iteration paths for each times-step with the solution, memory-consuming
        :param initial_guess_a_offset: a function of timestep (i), allows to choose an initial guess to find multiple solutions, see Remark 4.6. "None" means that the value will be taken from the previous time-step.
        :param iter_max: maximum number of iteration per time-step, stopping condition to detect divergence
        :param value_max: maximum norm of a solution, stopping condition to detect divergence
        :param eps: precision to detect convergence, stopping condition for the iterations
        :return: time-evolution of the state variables, a Solution object with many array fields
        """
        print("Starting the implicit catch-up method in V x R^m")
        T = np.zeros(nsteps + 1)
        YvA = np.zeros((self.dimV + self.m, nsteps + 1))
        Iters = np.zeros(nsteps)
        if save_iterations:
            iteration_steps=[]

        #seting up the description of the problem of the form VS @ [y_v, a] <= b_func(t,a), where [y_v, a] is the unknown
        #this particular description can be derived by substituting y = vbasis @ y_v into formula (66)
        self.box_offset_func = lambda t: self.G @ self.r(t) - self.F @ self.f(t)
        VS = np.block([[-self.vbasis, -self.S],
                       [self.vbasis, -self.S]])
        HS = self.H @ self.S
        b_func = lambda t, a: np.hstack((self.c0 - HS @ a - self.box_offset_func(t),
                                         self.c0 - HS @ a + self.box_offset_func(t)))

        self.moving_set_func = lambda t, a: Polytope(VS, b_func(t, a), None, None)

        self.S_V = np.block([[self.vbasis.T @ self.Kinv @ self.vbasis, np.zeros((self.dimV, self.m))],
                       [np.zeros((self.m, self.dimV)), np.eye(self.m)]])

        #initial condition
        yva0 = np.hstack((self.P_V_coords @ self.K @ e0 + self.G_V @ self.r(t0), a0))

        T[0] = t0
        YvA[:, 0] = yva0[:]
        quadprog_solver = McGibbonQuadprog() #choosing the optimization implementation for the projection step
        for i in range(nsteps):
            t_0 = T[i]
            yva_0 = YvA[:, i]
            t_1 = t_0 + dt

            if initial_guess_a_offset is None:
                yva_1 = np.copy(yva_0)
            else:
                yva_1 = yva_0 +  np.hstack((np.zeros(self.dimV),initial_guess_a_offset(i)))

            iter_counter = 0
            if save_iterations:
                iteration_steps.append(np.expand_dims(yva_1,axis=1))
            fixed_point_found = False

            while not fixed_point_found:
                C_1 = self.moving_set_func(t_1, yva_1[self.dimV:])
                try:
                    yva_1new = C_1.projection(self.S_V, yva_0, quadprog_solver)
                except ValueError:
                    raise NameError("Can't perform a step number " + str(i) + ", iteration " + str(iter_counter)+", projection failed")
                if np.linalg.norm(yva_1new, ord=np.inf) > value_max:
                    raise NameError("Can't perform a step number " + str(i) + ", iteration " + str(iter_counter)+", iterations diverged")

                sys.stdout.write(
                    "\r Completed iteration " + str(iter_counter) + ", step " + str(i + 1) + " of " + str(nsteps))
                sys.stdout.flush()

                fixed_point_found = (np.linalg.norm(yva_1new[self.dimV:]-yva_1[self.dimV:]) < eps)
                yva_1 = np.copy(yva_1new)

                iter_counter += 1
                if iter_counter > iter_max:
                    raise NameError("Can't perform a step number " + str(i) + ", maximum amount of iterations reached.")
                if save_iterations:
                    iteration_steps[i] = np.hstack((iteration_steps[i], np.expand_dims(yva_1, axis=1)))

            T[i + 1] = t_1
            YvA[:, i + 1] = yva_1
            Iters[i] = iter_counter

        print("")
        # end of implicit catch-up

        # produce all the variables of the solution
        Yv = np.zeros((self.dimV, nsteps + 1))
        Y = np.zeros((self.m, nsteps + 1))
        E = np.zeros((self.m, nsteps + 1))
        A = np.zeros((self.m, nsteps + 1))
        P = np.zeros((self.m, nsteps + 1))
        X = np.zeros((self.m, nsteps + 1))
        Sigma = np.zeros((self.m, nsteps + 1))  # stresses
        Zeta = np.zeros((self.n * self.d, nsteps + 1))  # displacements
        XI = np.zeros((self.n * self.d, nsteps + 1))  # positions
        Rho = np.zeros((self.q, nsteps + 1))  # reaction force values

        for i in range(nsteps + 1):  # can be further vectorized, provided that r and f are vectorized
            Yv[:, i] = YvA[0:self.dimV, i]
            Y[:, i] = self.vbasis @ Yv[:, i]
            Sigma[:, i] = Y[:, i] - self.box_offset_func(T[i])
            A[:, i] = YvA[self.dimV:, i]
            E[:, i] = self.Kinv @ Sigma[:, i]
            Rho[:, i] = self.Rp.T @ self.d_xi_phi.T @ (-Sigma[:, i] + self.F1 @ self.f(T[i]))
            if i == 0:
                P[:, i] = p0[:]
            else:
                P[:, i] = P[:, i - 1] + self.Sinv @ np.diag(np.sign(Sigma[:, i])) @ (
                        A[:, i] - A[:, i - 1])
            X[:, i] = E[:, i] + P[:, i]
            Zeta[:, i] = self.full_CM_pinv @ np.hstack((X[:, i], -self.r(T[i]) - self.R @ self.xi0))
            XI[:, i] = self.xi0 + Zeta[:, i]

        solution = self.Solution()
        solution.T = T
        solution.YvA = YvA
        solution.Y = Y
        solution.Yv = Yv
        solution.E = E
        solution.Sigma = Sigma
        solution.A = A
        solution.P = P
        solution.X = X
        solution.Zeta = Zeta
        solution.XI = XI
        solution.Rho = Rho
        solution.Iters = Iters
        if save_iterations:
            solution.iteration_steps=iteration_steps
        print("Finished the implicit catch-up method in V x R^m")
        return solution

    def solve_iteratively(self, e0, t0, dt, nsteps, a0, p0, save_iterations=False, initial_guess_a_offset=None,  iter_max = 5000, value_max = 1e+12, eps = 1e-12):
        """
        Implicit catch-up scheme, Algorithm 1
        :param e0: initial values for elastic elongation, m-vector
        :param t0: initial time
        :param dt: time-step
        :param nsteps: number of time-steps
        :param a0: initial values for the damage variable, m-vector
        :param p0: initial values for plastic elongation, m-vector
        :param save_iterations: set to True to include iteration paths for each times-step with the solution, memory-consuming
        :param initial_guess_a_offset: a function of timestep (i), allows to choose an initial guess to find multiple solutions, see Remark 4.6. "None" means that the value will be taken from the previous time-step.
        :param iter_max: maximum number of iteration per time-step, stopping condition to detect divergence
        :param value_max: maximum norm of a solution, stopping condition to detect divergence
        :param eps: precision to detect convergence, stopping condition for the iterations
        :return: time-evolution of the state variables, a Solution object with many array fields
        """
        print("Starting the implicit catch-up method in R^m x R^m")
        T = np.zeros(nsteps + 1)
        YA = np.zeros((2*self.m, nsteps + 1))
        Iters = np.zeros(nsteps)
        if save_iterations:
            iteration_steps=[]

        # seting up the description of the problem of the form IS @ [y, a] <= b_func(t,a), Aeq @ [y,a] = beq, where [y, a] is the unknown
        # this particular description can be derived from formula (66), Aeq follows from Proposition 4.1 i)
        self.box_offset_func = lambda t: self.G @ self.r(t) - self.F @ self.f(t)
        IS = np.block([[-np.eye(self.m), -self.S],
                       [np.eye(self.m), -self.S]])
        HS = self.H @ self.S
        b_func = lambda t, a: np.hstack((self.c0 - HS @ a - self.box_offset_func(t),
                                         self.c0 - HS @ a + self.box_offset_func(t)))

        Aeq = np.hstack(((self.Kinv @ self.ubasis).T, np.zeros((self.dimU, self.m))))
        beq = np.zeros(Aeq.shape[0])

        self.moving_set_func = lambda t, a: Polytope(IS, b_func(t, a), Aeq, beq)

        self.KinvI = np.block([[self.Kinv, np.zeros((self.m, self.m))],
                             [np.zeros((self.m, self.m)), np.eye(self.m)]])

        #initial condition
        ya0 = np.hstack((self.K @ e0 + self.box_offset_func(t0), a0))

        T[0] = t0
        YA[:, 0] = ya0[:]
        quadprog_solver = McGibbonQuadprog() #choosing the optimization implementation for the projection step
        for i in range(nsteps):
            t_0 = T[i]
            ya_0 = YA[:, i]
            t_1 = t_0 + dt

            if initial_guess_a_offset is None:
                ya_1 = np.copy(ya_0)
            else:
                ya_1 = ya_0 +  np.hstack((np.zeros(self.m),initial_guess_a_offset(i)))

            iter_counter = 0
            if save_iterations:
                iteration_steps.append(np.expand_dims(ya_1,axis=1))
            fixed_point_found = False

            while not fixed_point_found:
                C_1 = self.moving_set_func(t_1, ya_1[self.m:])
                try:
                    ya_1new = C_1.projection(self.KinvI, ya_0, quadprog_solver)
                except ValueError:
                    raise NameError("Can't perform a step number " + str(i) + ", iteration " + str(iter_counter)+", projection failed")
                if np.linalg.norm(ya_1new, ord=np.inf) > value_max:
                    raise NameError("Can't perform a step number " + str(i) + ", iteration " + str(iter_counter)+", iterations diverged")

                sys.stdout.write(
                    "\r Completed iteration " + str(iter_counter) + ", step " + str(i + 1) + " of " + str(nsteps))
                sys.stdout.flush()

                fixed_point_found = (np.linalg.norm(ya_1new[self.m:]-ya_1[self.m:]) < eps)
                ya_1 = np.copy(ya_1new)

                iter_counter += 1
                if iter_counter > iter_max:
                    raise NameError("Can't perform a step number " + str(i) + ", maximum amount of iterations reached.")
                if save_iterations:
                    iteration_steps[i] = np.hstack((iteration_steps[i], np.expand_dims(ya_1, axis=1)))

            T[i + 1] = t_1
            YA[:, i + 1] = ya_1
            Iters[i] = iter_counter

        print("")
        # end of implicit catch-up

        # produce all the variables of the solution
        Y = np.zeros((self.m, nsteps + 1))
        E = np.zeros((self.m, nsteps + 1))
        A = np.zeros((self.m, nsteps + 1))
        P = np.zeros((self.m, nsteps + 1))
        X = np.zeros((self.m, nsteps + 1))
        Sigma = np.zeros((self.m, nsteps + 1))  # stresses
        Zeta = np.zeros((self.n * self.d, nsteps + 1))  # displacements
        XI = np.zeros((self.n * self.d, nsteps + 1))  # positions
        Rho = np.zeros((self.q, nsteps + 1))  # reaction force values

        for i in range(nsteps + 1):  # can be further vectorized, provided that r and f are vectorized
            Y[:, i] = YA[0:self.m, i]
            Sigma[:, i] = Y[:, i] - self.box_offset_func(T[i])
            A[:, i] = YA[self.m:, i]
            E[:, i] = self.Kinv @ Sigma[:, i]
            Rho[:, i] = self.Rp.T @ self.d_xi_phi.T @ (-Sigma[:, i] + self.F1 @ self.f(T[i]))
            if i == 0:
                P[:, i] = p0[:]
            else:
                P[:, i] = P[:, i - 1] + self.Sinv @ np.diag(np.sign(Sigma[:, i])) @ (
                        A[:, i] - A[:, i - 1])
            X[:, i] = E[:, i] + P[:, i]
            Zeta[:, i] = self.full_CM_pinv @ np.hstack((X[:, i], -self.r(T[i]) - self.R @ self.xi0))
            XI[:, i] = self.xi0 + Zeta[:, i]

        solution = self.Solution()
        solution.T = T
        solution.YA = YA
        solution.Y = Y
        solution.E = E
        solution.Sigma = Sigma
        solution.A = A
        solution.P = P
        solution.X = X
        solution.Zeta = Zeta
        solution.XI = XI
        solution.Rho = Rho
        solution.Iters = Iters
        if save_iterations:
            solution.iteration_steps=iteration_steps
        print("Finished the implicit catch-up method in R^m x R^m")
        return solution