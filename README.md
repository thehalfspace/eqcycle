Spectral element method for earthquake cycle simulations with dynamic treatment of inertial effects.

development phase: working on better algorithms and parallelization.

Written in Julia 1.1

Adapted from Kaneko et al. (2011), and J.P. Ampuero's [SEMLAB](https://www.mathworks.com/matlabcentral/fileexchange/6154-semlab).

This project is incomplete, work in progress. To run the program 

1. change parameters in src/parameters/par.jl 
3. In the terminal, type `export JULIA_NUM_THREADS=4`. Alternatively, you can edit this in the scriptrun.jl line 11. When running on cluster, I write this in the pbs file and when running on computer, I edit in the script.
3. In line 23 of run.jl, P = setParameters(fzdepth, resolution), you can change the resolution as required. The target resolution is 24-32.
4. run the program in terminal using "julia run.jl"
