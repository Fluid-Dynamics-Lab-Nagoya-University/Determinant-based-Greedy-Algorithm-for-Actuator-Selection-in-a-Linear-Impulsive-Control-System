# Determinant-based-Greedy-Algorithm-for-Actuator-Selection-in-a-Linear-Impulsive-Control-System
This repository contains MATLAB R2024a codes for actuator placement optimization in a linear impulsive control system by determinant-based greedy algorithm.

## Codes
- P_act_opt_one_initial_cond.m: Main code to optimize actuator locations and compute the repsponse for one initial condition
- P_act_opt_multi_initial_cond.m: Main code to optimize actuator locations and compute the repsponses for multiple initial conditions
- F_act_dg_under1.m: Function for optimizing actuator locations
- F_act_dg_over1.m: Function for optimizing actuator locations
- F_norm_count2.m: Function for computing the 2 norm of the terminal output

## Note
It is needed to download "herdif.m" with relating functions from https://appliedmaths.sun.ac.za/~weideman/research/differ.html and put them inside a directory named DMSUITE to run the main codes.

## How to cite
Please cite the following paper.
M. Watanabe, Y. Sasaki, T. Nagata, K. Yamada, D. Tsubakino, J. Ito, and T. Nonomura, "Actuator placement optimization in a linear impulsive control system by determinant-based greedy algorithm", IEEE Access, Vol. 13, pp. 207581-207595, 2025, doi: 10.1109/ACCESS.2025.3641758. 
