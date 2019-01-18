Spectral element method for earthquake cycle simulations with fully dynamic inertial effects.

Adapted from Kaneko et al. (2011), and J.P. Ampuero's [SEMLAB](https://www.mathworks.com/matlabcentral/fileexchange/6154-semlab).

Running Instructions
====================

- open run.jl: 
    1.  Add the number of processors/workers in line 13. Default is 4.
    
    2. set resolution in line 30. Resolution is an integer. Default resolution = 1. Target resolution = 20. 
    3. Set the file name in line 45, run.jl.
- In the terminal run "julia run.jl". The output is a binary stored in data/ directory.

-----------------------------------------
Directory structure and File information
=========================================
eqcycle
    |- run.jl: The file for running the program
    |- output.jl: Data structure for reading the output files
    |- analyze_results.jl: De-serialize the binary output and get the variables for plotting.
    |- src: This folder contains the source files 
    |   |- gll_xwh: pre-tabulated interpolation points.
    |   |- paramters
    |   |   |-defaultParamters.jl: the default set of parameters.
    |   |- initialConditions
    |   |   |-defaultInitialConditions.jl: the default set of initial friction and stress conditions
    |   |- setup.jl: set up the mesh and initial properties. These are the variables not changing through the simulation.
    |    |- main.jl: The main function for the time-loop.
    |    |- Assemble.jl: Assemble the mass, stiffness, and other global matrices from the mesh.
    |    |- BoundaryMatrix.jl: set the boundary conditions
    |    |- dtevol.jl: compute the next timestep
    |    |- FindNearestNode.jl: Given an arbitrary point, find the closest point in the mesh.
    |    |- GetGLL.jl: read from pre-tabulated interpolation points.
    |    |- NRsearch.jl: Newton-Rhapson root finding algorithm.
    |    |- PCG.jl: preconditioned conjugate gradient method to invert the matrix.
    |    |- otherfunctions.jl: some other element wise computations.
    |- su|bmit_jl: Script for submitting jobs on the flux. Run on flux as "./submit_jl <jobname>"
    |- xfer_up: transfer output from flux to xfer server.
    |- xfer_down: transfer output from xfer server to my computer.
    |- scripts: scripts for plotting.


