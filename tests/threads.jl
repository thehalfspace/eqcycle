# Threaded matrix multiplication
import Base: eltype, size
#  import LinearAlgebra: A_mul_B!
using Base.Threads

struct ThreadedMul{Tv,Ti}
        A::SparseMatrixCSC{Tv,Ti}
end

function LinearAlgebra.mul!(y::AbstractVector, M::ThreadedMul, x::AbstractVector)
    @threads for i = 1 : M.A.n
         _threaded_mul!(y, M.A, x, i)
    end
     y
end

@inline function _threaded_mul!(y, A::SparseMatrixCSC{Tv}, x, i) where {Tv}
    s = zero(Tv)
    @inbounds for j = A.colptr[i] : A.colptr[i + 1] - 1
        s += A.nzval[j] * x[A.rowval[j]]
    end
    
    @inbounds y[i] = s
    y
end
eltype(M::ThreadedMul) = eltype(M.A)
size(M::ThreadedMul, I...) = size(M.A, I...)

    
    # P[1] = integer
    # P[2] = float
    # P[3] = float array
    # P[4] = integer array
    # P[5] = iglob
    # P[6] = ksparse
    # P[7] = H
    # P[8] = Ht
    # P[9] = damage_idx
    
    # damage evolution
    #  damage_idx = P[9]


    #  W_orig = W[:,:,damage_idx]
    #  damage_amount::Float64 = 1.0

    #  wgll2::Array{Float64,2} = S.wgll*S.wgll';
    
    # Time solver variables
    dt = P[2].dt0
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

    # dt modified slightly for damping
    if P[2].ETA != 0
        dt = dt/sqrt(1 + 2*P[2].ETA)
    end

    # Initialize kinematic field: global arrays
    d = zeros(P[1].nglob)
    v = zeros(P[1].nglob)
    v .= 0.5e-3
    a = zeros(P[1].nglob)
    
    #.....................................
    # Stresses and time related variables
    #.....................................
    tau = zeros(P[1].FltNglob)
    FaultC = zeros(P[1].FltNglob)
    Vf =  zeros(P[1].FltNglob)
    Vf1 = zeros(P[1].FltNglob)
    Vf2 = zeros(P[1].FltNglob)
    Vf0 = zeros(length(P[4].iFlt))
    FltVfree = zeros(length(P[4].iFlt))
    psi = zeros(P[1].FltNglob)
    psi0 = zeros(P[1].FltNglob)
    psi1 = zeros(P[1].FltNglob)
    psi2 = zeros(P[1].FltNglob)
    tau1 = zeros(P[1].FltNglob)
    tau2 = zeros(P[1].FltNglob)
    tau3 = zeros(P[1].FltNglob)
    tauAB = zeros(P[1].FltNglob)
    # Vectorized element calculations
    #  a_elem::Array{Float64,3} = zeros(size(P[5])) 
    #  Conn = sparse(P[5][:], collect(1:size(P[6][:])[1]),1)


    # Initial state variable
    psi = P[3].tauo./(P[3].Seff.*P[3].ccb) - P[3].fo./P[3].ccb - (P[3].cca./P[3].ccb).*log.(2*v[P[4].iFlt]./P[3].Vo)
    psi0 .= psi[:]

    isolver = 1
    
    # Skip lines 486-490
    # Skip lines 492-507: Outloc1, 2, variables.

    # Display important parameters
    #  println("Total number of nodes on fault: ", P.FltNglob)
    #  println("Average node spacing: ", P.LX/(P.FltNglob-1))
    #  @printf("dt: %1.09f s\n", dt)

    # Some more initializations
    r = zeros(P[1].nglob)
    beta_ = zeros(P[1].nglob)
    alpha_ = zeros(P[1].nglob)

    F = zeros(P[1].nglob)
    dPre = zeros(P[1].nglob)
    vPre = zeros(P[1].nglob)
    dd = zeros(P[1].nglob)
    dacum = zeros(P[1].nglob)
    dnew = zeros(length(P[4].FltNI))
    nseis = length(P[4].out_seis)

    # Preallocate variables with unknown size
    #  seismic_stress, seismic_slipvel, seismic_slip
    #  index_eq
    #  is_stress, is_slipvel, is_slip
    #  dSeis, vSeis, aSeis
    #  tStart, tEnd
    #  taubefore, tauafter, delfafter
    #  hypo, time_, Vfmax
    
    output = results(zeros(P[1].FltNglob, 150000), zeros(P[1].FltNglob, 150000), 
                     zeros(P[1].FltNglob, 150000), 
                     zeros(200000), 
                     zeros(P[1].FltNglob, 2000), zeros(P[1].FltNglob, 2000), 
                     zeros(P[1].FltNglob, 2000),
                     zeros(100000,nseis), zeros(100000,nseis), zeros(100000,nseis),
                     zeros(400), zeros(400), 
                     zeros(P[1].FltNglob, 400), zeros(P[1].FltNglob, 400), 
                     zeros(P[1].FltNglob, 400), zeros(400), zeros(500000), 
                     zeros(500000))
    
    tvsx = 2*P[1].yr2sec  # 2 years for interseismic period
    tvsxinc = tvsx

    tevneinc = 5    # 5 second for seismic period
    delfref = zeros(P[1].FltNglob)

    # Iterators
    idelevne= 0
    tevneb= 0
    tevne= 0
    ntvsx= 0
    nevne= 0
    slipstart= 0
    idd = 0
    it_s = 0; it_e = 0

    v = v[:] .- 0.5*P[2].Vpl
    Vf = 2*v[P[4].iFlt]
    iFBC = findall(abs.(P[3].FltX) .> 24e3)
    NFBC = length(iFBC) + 1
    Vf[iFBC] .= 0.


    v[P[4].FltIglobBC] .= 0.
    

    # on fault and off fault stiffness
    Ksparse = P[6]
    #  Ks2 = P[7]
    kni = -Ksparse[P[4].FltNI, P[4].FltNI]
    #  kni2 = -Ks2[P[4].FltNI, P[4].FltNI]

    nKsparse = -Ksparse

    # multigrid
    ml = ruge_stuben(kni)
    p = aspreconditioner(ml)
    tmp = copy(a)
    
    Ksparse = ThreadedMul(Ksparse)
    nKsparse = ThreadedMul(nKsparse)
    kni = ThreadedMul(kni)
    

    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0.

    #  while t < P[1].Total_time
        it = it + 1
        t = t + dt

        output.time_[it] = t 

        #  if isolver == 1

            vPre .= v
            dPre .= d

            Vf0 .= 2*v[P[4].iFlt] .+ P[2].Vpl
            Vf  .= Vf0

            #  for p1 = 1:2
                
                # Compute the on-Fault displacement
                F .= 0.
                F[P[4].iFlt] .= dPre[P[4].iFlt] .+ v[P[4].iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[P[4].FltNI]

                
                # Solve d = K^-1F by MGCG
                rhs = (mul!(tmp,Ksparse,F))[P[4].FltNI]
                
                # direct inversion
                #  dnew = -(kni\rhs)

                # mgcg
                dnew = cg!(dnew, kni, rhs, Pl=p, tol=1e-6)

                
                # update displacement on the medium
                d[P[4].FltNI] .= dnew

                # make d = F on the fault
                d[P[4].iFlt] .= F[P[4].iFlt]

                # Compute on-fault stress
                a .= 0.

                # Compute forcing (acceleration) for each element
                mul!(a,Ksparse,d)
                #  a = Ksparse*d

                tau1 .= -a[P[4].iFlt]./P[3].FltB
