# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).


import numpy as np

from solver.boundary_conditions import AffineBoundaryCondition

class Grid:
    """
    A class to handle 2d grid lattices
    """

    def __init__(self, n1, n2, is_node_func, xi0_func, add_springs_func, k_func, c0_func,
                 r_and_r_prime_component_functions, force_func, s_func=None, h_func=None):
        self.n1 = n1
        self.n2 = n2

        self.d = 2
        self.node_id_by_coord = -1*np.ones((n1, n2)).astype(int) #-1 means no node
        self.node_coord_by_id = []
        self.connections = []

        k = 0
        for i in range(n1):
            for j in range(n2):
                if is_node_func((i,j)):
                    self.node_id_by_coord[i, j] = k
                    self.node_coord_by_id.append((i,j))
                    k = k + 1
        self.n = k
        self.xi = np.zeros(self.n * self.d)
        self.Q = np.zeros((self.n,0))
        for k in range(self.n):
            (i, j) = self.node_coord_by_id[k]
            xi=xi0_func((i, j))
            self.xi[2*k] =  xi[0]
            self.xi[2*k+1] = xi[1]
            termins = add_springs_func((i, j))
            for k1 in range(len(termins)):
                (i1,j1)=termins[k1]
                edge_vect = np.zeros((self.n, 1))
                edge_vect[k,0] = 1
                edge_vect[self.node_id_by_coord[i1, j1], 0] = -1
                self.Q = np.append(self.Q, edge_vect, axis=1)
                self.connections.append((k,self.node_id_by_coord[i1,j1]))

        self.m = self.Q.shape[1]

        #setting up zero initial conditions
        self.e0 = np.zeros(self.m)
        self.a0 = np.zeros(self.m)
        self.p0 = np.zeros(self.m)
        self.zeta0 = np.zeros(self.n * self.d)

        #setting up the elastic and plastic parameters of the springs
        self.k=np.zeros(self.m)
        self.cminus = np.zeros(self.m)
        self.c0 = np.zeros(self.m)
        for i in range(self.m):
            self.k[i] = k_func(self.node_coord_by_id[self.connections[i][0]], self.node_coord_by_id[self.connections[i][1]])
            self.c0[i] = c0_func(self.node_coord_by_id[self.connections[i][0]], self.node_coord_by_id[self.connections[i][1]])

        if s_func is not None:
            self.s = np.zeros(self.m)
            for i in range(self.m):
                self.s[i] = s_func(self.node_coord_by_id[self.connections[i][0]],
                                   self.node_coord_by_id[self.connections[i][1]])
        if h_func is not None:
            self.h = np.zeros(self.m)
            for i in range(self.m):
                self.h[i] = h_func(self.node_coord_by_id[self.connections[i][0]],
                                   self.node_coord_by_id[self.connections[i][1]])

        self.boundary_condition = self.get_node_wise_boundary_condition(r_and_r_prime_component_functions, force_func)

    def get_node_wise_boundary_condition(self, r_and_r_prime_component_functions, force_func):
        """
        Combines the data on individual springs to return the data in the format of AffineBoundaryCondition
        :param r_and_r_prime_component_functions: the function defining the displacement boundary condition at a node
        :param force_func: the function defining an external force at a node
        :return: an instance of AffineBoundaryCondition for the entire lattice
        """
        q = 0
        R = np.zeros((0, self.n * self.d))
        self.rho_list = []
        for k in range(self.n):
            (i, j) = self.node_coord_by_id[k]
            for component in range(self.d):
                r_func = r_and_r_prime_component_functions((i,j), component) #checking that there is a constraint on that node
                if r_func is not None:
                #a new row of R which corresponds to the node and the component
                    R = np.vstack((R, np.zeros((1, self.n * self.d))))
                    R[q, k * self.d + component] = 1
                    self.rho_list.append(((i, j), component, r_func))
                    q = q + 1
        def r(t):
            r = np.zeros(q)
            counter = 0
            for constr in self.rho_list:
                r_func = constr[2]
                r[counter] = -r_func(t)[0]
                counter = counter + 1
            return r - R @ (self.xi + self.zeta0)

        def r_prime(t):
            r = np.zeros(q)
            counter = 0
            for constr in self.rho_list:
                r_func = constr[2]
                r[counter] = -r_func(t)[1]
                counter = counter + 1
            return r

        def f(t):
            forces = np.zeros(self.n * self.d)
            for k in range(self.n):
                f = force_func(self.node_coord_by_id[k])
                forces[self.d * k] = f(t)[0]
                forces[self.d * k + 1] = f(t)[1]
            return forces

        return AffineBoundaryCondition(q, R, r, r_prime, f)






