## Code to accompany: *[Most incompatible measurements for robust steering tests](https://arxiv.org/abs/1704.02994)*
#### Jessica Bavaresco, Marco Túlio Quintino, Leonardo Guerini, Thiago O. Maciel, Daniel Cavalcanti, and Marcelo Terra Cunha

This is a repository for all code which was written for the article "*Most incompatible measurements for robust steering tests*. J. Bavaresco, M. T. Quintino, L. Guerini, T. O. Maciel, D. Cavalcanti, and M. Terra Cunha. [*Phys. Rev. A* **96**, 022110 (2017),](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.022110) [arXiv:1704.02994 [quant-ph]](https://arxiv.org/abs/1704.02994)."

All code is written in MATLAB and requires:
- [Yalmip](https://yalmip.github.io) - a free MATLAB toolbox for rapid prototyping of optimization problems
- [MOSEK](https://www.mosek.com) - a software package for solving mathematical optimization problems (under the free personal academic license)
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory

The code is separeted in three main methods that perform different tasks: the search algorithm, the see-saw algorithm, and the outer polytope approximation.

- Search algorithm: 
The search algorithm is employed in optimizations over sets of measurements that are restricted in the number of measurements, number of outcomes of each measurement, and also possibly in the structure of the POVMs. The code made available here is specifically for the case of qubit measurements.

  - [search_qubit_planproj](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/search_qubit_planproj.m): calculates upper bounds for the critical visibility of two-qubit quantum states subjected to N planar projective measurements and provides a candidate for the optimal set of measurements.
  
  - [search_qubit_proj](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/search_qubit_proj.m): calculates upper bounds for the critical visibility of two-qubit quantum states subjected to N projective measurements and provides a candidate for the optimal set of measurements.
  
  - [search_qubit_trine](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/search_qubit_trine.m): calculates upper bounds for the critical visibility of two-qubit quantum states subjected to N regular trine measurements and provides a candidate for the optimal set of measurements.
  
  - [search_qubit_tetra](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/search_qubit_tetra.m): calculates upper bounds for the critical visibility of two-qubit quantum states subjected to N regular tetrahedron measurements and provides a candidate for the optimal set of measurements.
  
  - [search_qubit_genpovm](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/search_qubit_genpovm.m): calculates upper bounds for the critical visibility of two-qubit quantum states subjected to general N k-outcome measurements and provides a candidate for the optimal set of measurements.
  
- See-saw algorithm:
The see-saw algorithm optimizes over sets of general POVMs that are restricted only in the number of measurements and in the number of outcomes for each measurement. The code made available here can be used for states of any dimension.

  - [seesaw_wnr](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/seesaw_wnr.m): see-saw method for calculating an upper bound on the critical visibility to white noise of a given quantum state subjected to N k-outcome general POVMs. Provides a candidate to the optimal set of measurements to steer the given quantum state in the specified scenario.
  
  - [seesaw_gr](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/seesaw_gr.m): see-saw method for calculating a lower bound on generalized robustness of steering of a given quantum state subjected to N k-outcome general POVMs. Provides a candidate to the optimal set of measurements to steer the given quantum state in the specified scenario.
 
- Outer polytope approximation<sup>§</sup>:
The outer polytope approximation method estimates the amount of white noise that must me applied to a quantum state in order to guarantee that it is unsteerable under the specified scenario. The code made available here is specifically for the case of N=5 qubit projective measurements.

  - [polyapprox_qubit_5planproj](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/polyapprox_qubit_planproj.m): calculates lower bounds for the critical visibility of two-qubit quantum states subjected to N=5 (can be simply modified to any other number) planar projective measurements and provides LHS models for the assemblages in the set.
  
  - [polyapprox_qubit_5proj](https://github.com/jessicabavaresco/most-incompatible-measurements/blob/master/polyapprox_qubit_proj.m): calculates lower bounds for the critical visibility of two-qubit quantum states subjected to N=5 (can be simply modified to any other number) projective measurements and provides LHS models for the assemblages in the set.
 
<sup>§</sup>: The outer polytope approximation code also requires:
- [CDDMEX](http://control.ee.ethz.ch/~cohysys/cdd.php) - a free MATLAB interface for the CDD solver, a software package of efficient algorithms for polytope manipulation
