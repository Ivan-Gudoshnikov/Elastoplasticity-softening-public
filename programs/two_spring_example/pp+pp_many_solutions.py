# (C) Ivan Gudoshnikov, 2025
# See "Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"
# This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

from example0_projection_chart import example0_projection_chart
from solver.elastoplastic_process_linearized import ElastoplasticProcessLinearized
import numpy as np
import matplotlib.pyplot as plt

d = 1

Q = np.array([[1, 0],
              [-1, 1],
              [0, -1]])

xi0 = np.array([0., 1., 2.])
n = Q.shape[0]
m = Q.shape[1]

k = np.ones(m)
s = np.array([1.,1.])
h = np.array([1., 1.])

c0 = 0.01 * np.array([1., 1.])

q = 2

R = np.array([[1., 0., 0.],
              [0., 0., 1.]])

r_prime_vect = np.array([0, -1.])
r_prime = lambda t: r_prime_vect
r = lambda t: r_prime_vect * t - R @ xi0
f = lambda t: np.zeros(3)

t0 = 0
dt = 0.0004
nsteps = 10
process = ElastoplasticProcessLinearized(Q, xi0, k, c0, s, h, d, R, r, f, True)

t_ref = 0
e0 = c0
a0 = np.zeros(m)
p0 = np.zeros(m)

def initial_guess_offet(i):
    #return np.array([8.e-5, -1.e-5])  #try this to get another solution
    return np.array([0., 0.])

solution = process.solve_iteratively_in_V(e0, t0, dt, nsteps, a0, p0, save_iterations=True, initial_guess_a_offset=initial_guess_offet)

figA, axA = plt.subplots()
axA.plot(solution.T, solution.A[:, :].T, label=['a_1', 'a_2'])
axA.legend()
axA.set(title='Damage variable a')
plt.xlabel('t')
plt.ylabel('a')

i_to_examine = 0
_ = example0_projection_chart(process, solution, i_to_examine, 7e-4, 7e-4, draw_iterations=True)

plt.show()
