Spectral element method for earthquake cycle simulations with dynamic treatment of inertial effects.

Changelog:
    - Implemented [algebraic multigrid](https://github.com/JuliaLinearAlgebra/AlgebraicMultigrid.jl) conjugate gradient for the quasi-static solver.
    - Implemented multithreading parallelism

development phase: working on better algorithms and parallelization.

Written in Julia 1.0

Adapted from Kaneko et al. (2011), and J.P. Ampuero's [SEMLAB](https://github.com/jpampuero/semlab)

This project is incomplete, work in progress. To run the program 

1. change parameters in par.jl 
2. edit the initial conditions in src/initialConditions/defaultInitialConditions
3. Set the resolution and output files in run.jl
4. run the program in terminal using "julia run.jl"
5. Wait for eternity before the results pop up! (I am still working on making it faster)
