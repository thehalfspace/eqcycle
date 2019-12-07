#################################
# Run the simulations from here
#################################

# 1. Go to par.jl and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder

#  `export JULIA_NUM_THREADS=1`

using Printf
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using AlgebraicMultigrid
#  using BenchmarkTools
using StaticArrays
using IterativeSolvers
using FEMSparse
#  using CuthillMcKee
#  using Base.Threads
#  BLAS.set_num_threads(6)

include("$(@__DIR__)/par.jl")	    #	Set Parameters

P = setParameters(24e3,20)      # args = fault zone depth, resolution

include("$(@__DIR__)/src/dtevol.jl")          
include("$(@__DIR__)/src/NRsearch_serial.jl")
include("$(@__DIR__)/src/otherFunctions_serial.jl")

include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed O = @time main(P)

#  description = "homogeneous medium with high resolution"

# Save output to file
using Serialization
open("$(@__DIR__)/data_new/width75_cs60.out", "w") do f
    serialize(f,O)
    serialize(f, simulation_time)
    serialize(f, P)
end

println("\n")

@info("Simulation Complete!");
