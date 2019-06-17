# Distribute over workers
include("../src/damageEvol.jl")

function setParameters(FZdepth, res)

    LX::Int = 48e3  # depth dimension of rectangular domain
    LY::Int = 30e3 # off fault dimenstion of rectangular domain

    NelX::Int = 30*res # no. of elements in x
    NelY::Int = 20*res # no. of elements in y

    dxe::Float64 = LX/NelX #	Size of one element along X
    dye::Float64 = LY/NelY #	Size of one element along Y
    Nel::Int = NelX*NelY # Total no. of elements
    
    P::Int = 4		#	Lagrange polynomial degree
    NGLL::Int = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob::Int = NelX*(NGLL - 1) + 1

    # Jacobian for global -> local coordinate conversion
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    coefint1::Float64 = jac/dx_dxi^2
    coefint2::Float64 = jac/dy_deta^2

    #..................
    # TIME PARAMETERS
    #..................

    yr2sec::Int = 365*24*60*60
    
    Total_time::Int = 200*yr2sec     # Set the total time for simulation here

    CFL::Float64 = 0.6	#	Courant stability number
     
    IDstate::Int = 2    #   State variable equation type

    # Some other time variables used in the loop
    dtincf::Float64 = 1.2
    gamma_::Float64 = pi/4
    dtmax::Int = 50 * 24 * 60*60		# 5 days


    #...................
    # MEDIUM PROPERTIES
    #...................
    
    # default
    rho1::Float64 = 2670
    vs1::Float64 = 3464

    # The entire medium has low rigidity
    #  rho1::Float64 = 2500
    #  vs1::Float64 = 0.6*3464

    rho2::Float64 = 2500
    vs2::Float64 = 0.6*vs1
    
    ETA = 0.

    # Low velocity layer dimensions
    ThickX::Float64 = LX - ceil(FZdepth/dxe)*dxe # ~FZdepth m deep
    ThickY::Float64 = ceil(0.25e3/dye)*dye   # ~ 0.25*2 km wide

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading

    fo::Vector{Float64} 	= repeat([0.6], FltNglob)		#	Reference friction coefficient
    Vo::Vector{Float64} 	= repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'
    xLf::Vector{Float64} = repeat([0.008], FltNglob)#	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.001
    Vevne::Float64 = Vthres


    #-----------#
    #-----------#
    # SETUP
    #-----------#
    #-----------#
    
    #....................
    # 2D Mesh generation
    #....................
    iglob::Array{Int,3}, x::Vector{Float64}, y::Vector{Float64} = 
                        MeshBox!(NGLL, Nel, NelX, NelY, FltNglob, dxe, dye)
    x = x .- LX
    nglob::Int = length(x)

    # The derivatives of the Lagrange Polynomials were pre-tabulated 
    # xgll = location of the GLL nodes inside the reference segment [-1,1]
    xgll::Vector{Float64}, wgll::Vector{Float64}, H::Matrix{Float64} = GetGLL(NGLL)
    Ht::Matrix{Float64} = H'
    wgll2::Matrix{Float64} = wgll*wgll'

    #.................
    # Initialization
    #.................

    # For internal forces
    W::Array{Float64,3} = zeros(NGLL, NGLL, Nel)

    # Global Mass Matrix
    M::Vector{Float64} = zeros(nglob)

    # Mass+Damping matrix
    MC::Vector{Float64} = zeros(nglob)

    # Assemble mass and stiffness matrix
    M, W, dt::Float64, muMax = assemble!(NGLL, NelX, NelY, dxe, dye, ThickX,
                                ThickY, rho1, vs1, rho2, vs2, iglob, 
                                M, W, x, y, jac)

    return damageEvol!(NGLL, NelX, NelY, dxe, dye, ThickX,
                                ThickY, rho1, vs1, rho2, vs2, iglob, 
                                M, W, x, y, jac) 
    
    # Time solver variables
    dt = CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2
    
    #......................
    # Boundary conditions :
    #......................

    # Left boundary
    BcLC::Vector{Float64}, iBcL::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'L')

    # Right Boundary = free surface: nothing to do
    #  BcRC, iBcR = BoundaryMatrix(P, wgll, iglob, 'R')

    # Top Boundary
    BcTC::Vector{Float64}, iBcT::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'T')

    # Mass matrix at boundaries
    #  Mq = M[:]
    M[iBcL] .= M[iBcL] .+ half_dt*BcLC
    M[iBcT] .= M[iBcT] .+ half_dt*BcTC
    #  M[iBcR] .= M[iBcR] .+ half_dt*BcRC


    # Dynamic fault at bottom boundary
    FltB::Vector{Float64}, iFlt::Vector{Int} = BoundaryMatrix!(NGLL, NelX, NelY, rho1, vs1, rho2, vs2, dy_deta, dx_dxi, wgll, iglob, 'B')

    FltZ::Vector{Float64} = M[iFlt]./FltB /half_dt * 0.5
    FltX::Vector{Float64} = x[iFlt]
    
    #......................
    # Initial Conditions
    #......................
    cca::Vector{Float64}, ccb::Vector{Float64} = fricDepth(FltX)   # rate-state friction parameters
    Seff::Vector{Float64} = SeffDepth(FltX)       # effective normal stress
    tauo::Vector{Float64} = tauDepth(FltX)        # initial shear stress

    # Kelvin-Voigt Viscosity
    Nel_ETA::Int = 0
    if ETA !=0
        Nel_ETA = NelX
        x1 = 0.5*(1 .+ xgll')
        eta_taper = exp.(-pi*x1.^2)
        eta = ETA*dt*repeat([eta_taper], NGLL)

    else
        Nel_ETA = 0
    end

    # Compute XiLF used in timestep calculation
    XiLf::Vector{Float64} = XiLfFunc!(LX, FltNglob, gamma_, xLf, muMax, cca, ccb, Seff)
 
    # Find nodes that do not belong to the fault
    FltNI::Vector{Int} = deleteat!(collect(1:nglob), iFlt)
    
    # Compute diagonal of K
    diagKnew::Vector{Float64} = KdiagFunc!(FltNglob, NelY, NGLL, Nel, coefint1, coefint2, iglob, W, H, Ht, FltNI)
    
    # Fault boundary: indices where fault within 24 km
    fbc = reshape(iglob[:,1,:], length(iglob[:,1,:])) 
    idx = findall(fbc .== findall(x .== -24e3)[1] - 1)[1]    
    FltIglobBC::Vector{Int} = fbc[1:idx]

    # Display important parameters
    println("Total number of nodes on fault: ", FltNglob)
    println("Average node spacing: ", LX/(FltNglob-1))
    @printf("dt: %1.09f s\n", dt)

end

