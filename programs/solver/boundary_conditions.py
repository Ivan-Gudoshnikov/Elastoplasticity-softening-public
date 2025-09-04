# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

import numpy as np

class AffineBoundaryCondition:
    """
    A general class to hold the data on:
    - affine displacement boundary R(xi_0+ zeta) +r(t) == 0, see (LSM5) in the paper
    - an external force
    """
    def __init__(self, q, R, r, r_prime, f):
        self.q=q
        self.R=R
        self.r=r
        self.r_prime=r_prime
        self.f=f