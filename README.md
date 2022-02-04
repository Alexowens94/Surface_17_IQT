# Surface_17
Surface 17 sims for ion trapping quantum computing architecture 

Requires ProjectQ (installed with it's fast C++ circuit simulator) 
(get Visual C++ Build Tools and the Microsoft Windows SDK prior to doing installing projectQ via pip)

*** This project has grown horrendously untidy - we are in the process of refactoring an object oriented version to facilitate the addition of complex features *** 

This repo is for simulating quantum error correction, specifically the performance of the 17 qubit surface code under a variety of experimentally motivated error models. 
The idea is to use the logical error rates, pseudo-thresholds and information about the most damaging physical errors produced by these simulations to assess the suitability 
of design choices for use in fault-tolerant architectures. 
For example, in ion-traps the control errors on the 2-qubit gates depend on the occupation of the motional mode used to mediate the gate interaction - the motional
modes heat up during an algorithm and if they aren't cooled intermitently the gate error rate will also increase. Cooling operations take much longer than gate 
operations, and recooling frequently will result in much longer algorithm runtimes (see https://arxiv.org/pdf/2003.01293.pdf), 
providing more opportunity for idle qubits to dephase. Simulating both scenarios under realistic environmental noise levels can indicate the regimes in which cooling
should be advantageous or detrimental to circuit performance. 

Example outputs:
[Logical_qubit_decay](https://github.com/QECsims/Surface_17/blob/master/img/PDD_deph=0_p=0.01.png)
[Error_subsets](https://github.com/QECsims/Surface_17/blob/master/img/contributions_to_total_error_xyxx.png)
[Pl_no_cooling](https://github.com/QECsims/Surface_17/blob/master/img/pl_PDD_linear_heat.png)
[Pl_cooling](https://github.com/QECsims/Surface_17/blob/master/img/pl_cool.png)
[Dressed_qubit_vary_dephase_rate](https://github.com/QECsims/Surface_17/blob/master/img/dress_leaked_reg_with_reset_every_run.png)

The script importance_sampling_syndrome_processing.py provides an example running through 
1) calculating significant error subsets requiring simulation
2) simulating logical error rates given an error subset by randomly placing fixed number of errors according to error 
   model

plotting_imp_samp.py combines subset weight data with logical error rates to produce plots of the total error rate 
as in the above examples

The small circuit size allows for simulating 'advanced' error models including correlated and coherent errors, error rates varying in time/space, as well as 
leakage/loss from the code space. So far the models have been extended to include different dynamical decoupling schemes, different means of addressing leakage, 
circuit compilations in terms of different native gates and sympathetic cooling.
Simple Pauli models can easily be added as entries in the error_model.py file following the syntax there.

ProjectQ is used for the circuit simulations, and the methods for simulating noise are based on the work in: https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf 
The importance sampling method is explained further in: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.032341

The decoder employed uses a lookup table generated in generate_lookup_table.py and minimal syndrome processing rules
as described in https://iopscience.iop.org/article/10.1088/1367-2630/aab341/pdf (decoder is modular, and can be replaced with a differently weighted variant)

To discuss the repo reach me on A.Owens@sussex.ac.uk 
I don't intend to make the noise modelling underpinning the error models public, but happy to discuss via email.
https://arxiv.org/abs/2111.01913 is a good resource for a general laser-free *gate based* error model
https://journals-aps-org.ezproxy.sussex.ac.uk/prl/pdf/10.1103/PhysRevLett.111.140501 and https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.117.220501
outline the 'dressed state' dynamically decoupled qubit prominent in the simulationss along with some consideration of error sources. 

*********Importance Sampling*********
To make use of importance sampling with the chosen syndrome processing rules requires some ammendment to the protocol to account for the fact the circuit length
isn't known a priori. 
We have to simulate the error rate GIVEN the code would stop after a single syndrome measurement and the error rate given it would run for 2 syndrome measurements,
to determine the full logical error rate as below:
P_l = P(logic error | (n_errors AND s1))*P( (n_errors AND s1) )  + P(logic error | (n_errors AND s2))*P( (n_errors AND s2) ) 
where n_errors is a list of indicies indicating the total number of errors of each type occuring in the simulation, and s1/s2 indicate if the run continues for 
1 or 2 syndrome measurements. It is worth noting P( (n_errors AND s2) )  requires data about the probability of reaching s2 given a fixed number of errors in s1, 
which we determine during the simulation of s1 (i.e. P( (n_errors AND s2) ) can't be calculated purely combinatorically like P( (n_errors AND s1) ) can,
and requires some information from the simulations)

