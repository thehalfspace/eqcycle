include("$(@__DIR__)/src/GetGLL.jl")		 #	Polynomial interpolation
include("$(@__DIR__)/src/MeshBox.jl")		 # 	Build 2D mesh
include("$(@__DIR__)/src/Assemble.jl")       #   Assemble mass and stiffness matrix
include("$(@__DIR__)/src/Kassemble.jl")      #   Assemble mass and stiffness matrix
#  include("$(@__DIR__)/trapezoidFZ/Assemble.jl") #   Gaussian fault zone assemble
include("$(@__DIR__)/src/BoundaryMatrix.jl")    #	Boundary matrices
include("$(@__DIR__)/src/FindNearestNode.jl")   #	Nearest node for output
include("$(@__DIR__)/src/initialConditions/defaultInitialConditions.jl")


using Printf
using LinearAlgebra
using DelimitedFiles
using SparseArrays
using AlgebraicMultigrid
using BenchmarkTools
using StaticArrays
using IterativeSolvers




    LX = 48e3  # depth dimension of rectangular domain
    LY = 30e3 # off fault dimenstion of rectangular domain

    NelX = 30*res # no. of elements in x
    NelY = 20*res # no. of elements in y

    dxe = LX/NelX #	Size of one element along X
    dye = LY/NelY #	Size of one element along Y
    Nel = NelX*NelY # Total no. of elements
    
    P = 4		#	Lagrange polynomial degree
    NGLL = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob = NelX*(NGLL - 1) + 1

    # Jacobian for global -> local coordinate conversion
    dx_dxi = 0.5*dxe
    dy_deta = 0.5*dye
    jac = dx_dxi*dy_deta
    coefint1 = jac/dx_dxi^2
    coefint2 = jac/dy_deta^2

    #..................
    # TIME PARAMETERS
    #..................

    yr2sec = 365*24*60*60
    
    Total_time = 100*yr2sec     # Set the total time for simulation here

    CFL = 0.6	#	Courant stability number
     
    IDstate = 2    #   State variable equation type

    # Some other time variables used in the loop
    dtincf = 1.2
    gamma_ = pi/4
    dtmax = 50 * 24 * 60*60		# 5 days


    #...................
    # MEDIUM PROPERTIES
    #...................
    
    # default
    rho1 = 2670
    vs1 = 3464

    # The entire medium has low rigidity
    #  rho1 = 2500
    #  vs1 = 0.6*3464

    rho2 = 2500
    vs2 = 0.6*vs1
    
    ETA = 0.

    # Low velocity layer dimensions
    ThickX = LX - ceil(FZdepth/dxe)*dxe # ~FZdepth m deep
    ThickY = ceil(0.25e3/dye)*dye   # ~ 0.25*2 km wide

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl = 35e-3/yr2sec	#	Plate loading

    fo = repeat([0.6], FltNglob) #	Reference friction coefficient
    Vo = repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'
    xLf = repeat([0.008], FltNglob)    #	Dc (Lc) = 8 mm

    Vthres = 0.001
    Vevne = Vthres

    #-----------#
    #-----------#
    # SETUP
    #-----------#
    #-----------#
    
    #....................
    # 2D Mesh generation
    #....................
    iglob, x, y = MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)
    x = x .- LX
    nglob = length(x)

    # The derivatives of the Lagrange Polynomials were pre-tabulated 
    # xgll = location of the GLL nodes inside the reference segment [-1,1]
    xgll, wgll, H = GetGLL(NGLL)
    Ht = H'
    wgll2 = wgll*wgll'

    #.............................
    #   OUTPUT RECEIVER LOCATIONS
    #.............................
    # For now, it saves slip, sliprate, and stress at the nearest node specified.
    # My coordinates are weird, might change them later.
    # x coordinate = along dip fault length (always -ve below the free surface)
    # y coordinate = off-fault distance (+ve)
     

    x_out = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0].*(-1e3)  # x coordinate of receiver
    y_out = [66.0, 130.0, 198.0, 250.0, 330.0, 396.0]     # y coordinate of receiver
    #  n_receiver = length(x_receiver) # number of receivers

    x_out, y_out, out_seis, dist = FindNearestNode(x_out, y_out, x, y) 
    
    
    #.................
    # Initialization
    #.................

    # For internal forces
    #  W::Array{Float64,3} = zeros(NGLL, NGLL, Nel)

    # Global Mass Matrix
    M = zeros(nglob)

    # Mass+Damping matrix
    #  MC = zeros(nglob)

    # Assemble mass and stiffness matrix
    M, dt, muMax, damage_idx = assemble!(NGLL, NelX, NelY, dxe, dye, 
                        ThickX,ThickY, rho1, vs1, rho2, vs2, iglob,M, x, y, jac)

    # Stiffness Assembly
    #  Ksparse = stiffness_assembly(NGLL, NelX, NelY, 
                    #  nglob, dxe, dye, ThickX, ThickY, rho1, vs1, rho2, vs2, iglob)
