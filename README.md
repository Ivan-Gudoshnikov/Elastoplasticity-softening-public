# Elastoplasticity-softening-public

This is a Supplemental Material for the paper 
*"Elastoplasticity with softening as a state-dependent sweeping process: non-uniqueness of solutions and emergence of shear bands in lattices of springs"* 
by I. Gudoshnikov (Institute of Mathematics of the Czech Academy of Sciences, Žitná 25, 110 00, Praha 1, Czech Republic).

Included are full videos of the simulations of the Lattice Spring Models, obtained from the state-dependent sweeping process formulation via the implicit catch-up method (Algorithm 2 in the paper). 
Included are the examples from Section 6, as well as a few additional relevant examples.

Note: in the Figs. 10-16 a-h in the paper the color was adjusted so that red color means the maximum value of the respective variable on time-interval [0, 6]. In the supplemental videos, showing full evolution on time-interval [0,10], red color indicates the maximum value on [0,10]. 

The computer codes of the simulations are also included.

This research is supported by the Czech Science Foundation project GA24-10586S and the Czech Academy of Sciences (RVO: 67985840).

**A symmetric rectangular lattice with softening**
See Section 6.1.1 and Fig. 10 in the paper.

* square_10x15_softening_damage.mp4  - color indicates the evolution of the damage variable (a), thickness indicates a yielding spring (ongoing plastic  deformation).
* square_10x15_softening_stress.mp4  - color indicates the evolution of the stress variable (\sigma), thickness indicates a yielding spring (ongoing plastic deformation).


**A symmetric rectangular lattice with softening - the second solution**
See Section 6.1.1 and Fig. 3 in the paper.

* square_10x15_softening_initial_guess2_damage.mp4
* square_10x15_softening_initial_guess2_stress.mp4
* square_10x15_softening_initial_guess2_iters.png  - the graph of the number of iteration per time-step
* square_10x15_softening_initial_guess2_total_stress.png  - the graph of the total stress in the lattice (formula (112)) 


**A rectangular lattice with a defect with softening**
See Section 6.1.2 and Fig. 11 in the paper.

* square_10x15_softening_defect_to_right_damage.mp4
* square_10x15_softening_defect_to_right_stress.mp4


**A symmetric rectangular lattice with perfect plasticity**
See Section 6.1.3 and Fig. 12 in the paper.
  
* square_10x15_pp_damage.mp4
* square_10x15_pp_stress.mp4


**A rectangular lattice with a defect with perfect plasticity**
See Section 6.1.4 and Fig. 13 in the paper.

* square_10x15_pp_defect_to_right_damage.mp4
* square_10x15_pp_defect_to_right_stress.mp4


**A symmetric rectangular lattice with hardening**
See Section 6.1.5 and Fig. 14 in the paper.
  
* square_10x15_hardening_damage.mp4
* square_10x15_hardening_stress.mp4


**A rectangular lattice with a defect with hardening**
See Section 6.1.6 and Fig. 15 in the paper.

* square_10x15_hardening_defect_to_right_damage.mp4
* square_10x15_hardening_defect_to_right_stress.mp4


**A symmetric triangular lattice with softening**
Parameters of the springs as in Section 6.2 of the paper, the results are not included in the text.

* triangles_softening_damage.mp4
* triangles_softening_stress.mp4
* triangles_softening_iters.png
* triangles_softening_total_stress.png

**A triangular lattice with a defect with softening**
See Section 6.2 and Fig. 16 in the paper.

* triangles_softening_defect_damage.mp4
* triangles_softening_defect_stress.mp4

