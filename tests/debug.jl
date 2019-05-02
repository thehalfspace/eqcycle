function element_computation!(P::params_float, iglob::Array{Int,3}, F_local::Array{Float64}, H::Array{Float64,2}, Ht::Array{Float64,2}, W::Array{Float64,3}, Nel)
    a_elem = zeros(size(iglob))
        for eo in 1:Nel
            ig = iglob[:,:,eo]
            Wlocal = W[:,:,eo]
            locall = F_local[ig]
            a_elem[:,:,eo] =  P.coefint1*H*(Wlocal.*(Ht*locall)) + P.coefint2*(Wlocal.*(locall*H))*Ht
        end
    return a_elem
end

# Assemble K matrix as local parts [Nel,ngll,ngll]
function Kassemble(Nel, iglob, P, W, H, Ht)
    a_elem = zeros(size(iglob))
        for eo in 1:Nel
            Wlocal = W[:,:,eo]
            a_elem[:,:,eo] =  P.coefint1*H*(Wlocal.*Ht) + P.coefint2*(Wlocal.*H)*Ht
        end
    return a_elem
end

function kass2(kloc, iglob,nglob, nel)
    kglob = zeros(nglob, nglob)

    for e in 1:nel
        ig = iglob[:,:,e]
        kglob[ig[:],ig[:]] += repeat(kloc[ig][:],1,9)
    end

    kglob
end

#a_elem[:,:,eo] =  P.coefint1*H*(Wlocal.*(Ht*locall)) + P.coefint2*(Wlocal.*(locall*H))*Ht

W = P[6]
P2 = P[2];
iglob = P[5];
H = P[7]
Ht = P[8]
Nel = P[1].Nel
nglob = P[1].nglob

ngll = 5

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
    a_elem = zeros(size(P[5])) 
    Conn = sparse(P[5][:], collect(1:size(P[6][:])[1]),1)


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

    # Preallocate variables with unknown size
    output = results(zeros(P[1].FltNglob, 800000), zeros(P[1].FltNglob, 800000), 
                     zeros(P[1].FltNglob, 800000), zeros(1000000), 
                     zeros(P[1].FltNglob, 10000), zeros(P[1].FltNglob, 10000), 
                     zeros(P[1].FltNglob, 10000),zeros(1000), zeros(10000), 
                     zeros(P[1].FltNglob, 10000), zeros(P[1].FltNglob, 10000), 
                     zeros(P[1].FltNglob, 10000), zeros(10000), zeros(10000000), 
                     zeros(10000000))
    
    # Save output variables at certain timesteps: define those timesteps
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

            @inbounds @fastmath Vf0 .= 2*v[P[4].iFlt] .+ P[2].Vpl
            Vf  .= Vf0

            #  @inbounds for p1 = 1:2
                
                # Compute the on-Fault displacement
                F .= 0.
                @inbounds @fastmath F[P[4].iFlt] .= dPre[P[4].iFlt] .+ v[P[4].iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[P[4].FltNI]
                
                # Solve d = K^-1F by PCG
                dnew = PCG!(P[2], P[1].Nel, P[3].diagKnew, dnew, F, P[4].iFlt, 
                            P[4].FltNI,P[7], P[8], P[5], P[1].nglob, P[6], 
                            a_elem, Conn)#, a_local, dd_local, p_local)
                
                # update displacement on the medium
                d[P[4].FltNI] .= dnew

                #  # make d = F on the fault
                @inbounds @fastmath d[P[4].iFlt] .= F[P[4].iFlt]

                #  # Compute on-fault stress
                a .= 0.

                #  # Compute forcing (acceleration) for each element

                #  Threads.@threads for eo in 1:P[1].Nel
                    #  ig = P[5][:,:,eo]
                    #  a_elem[:,:,eo] = element_computation!(P[2], ig, d[ig], 
                                                          #  P[7], P[8], P[6][:,:,eo]) 
                #  end
                #  a= Conn*a_elem[:]

                #  a[P[4].FltIglobBC] .= 0.
                #  @inbounds @fastmath tau1 .= -a[P[4].iFlt]./P[3].FltB
                
                #  # Function to calculate sliprate
                #  Threads.@threads for j = NFBC:P[1].FltNglob
                    #  psi1[j], Vf1[j] = slrFunc!(P[3].xLf[j], P[3].Vo[j], P[3].fo[j], psi[j], psi1[j], Vf[j], Vf1[j], P[1].IDstate, tau1[j], P[3].tauo[j], tauAB[j], P[3].Seff[j], P[3].cca[j], P[3].ccb[j], dt)
                #  end

                #  Vf1[iFBC] .= P[2].Vpl
                #  @inbounds @fastmath Vf .= (Vf0 + Vf1)/2
                #  @inbounds @fastmath v[P[4].iFlt] .= 0.5*(Vf .- P[2].Vpl)

            #  end

            #  psi .= psi1[:]
            #  tau .= tau1[:]
            #  tau[iFBC] .= 0.
            #  Vf1[iFBC] .= P[2].Vpl

            #  @inbounds @fastmath v[P[4].iFlt] .= 0.5*(Vf1 .- P[2].Vpl)
            #  @inbounds @fastmath v[P[4].FltNI] .= (d[P[4].FltNI] .- dPre[P[4].FltNI])/dt

            #  # Line 731: P_MA: Omitted
            #  a .= 0.
            #  d[P[4].FltIglobBC] .= 0.
            #  v[P[4].FltIglobBC] .= 0.
            
            #  # If isolver != 1, or max slip rate is < 10^-2 m/s
        #  else
            
            #  dPre .= d
            #  vPre .= v

            #  # Update
            #  @inbounds @fastmath d .= d .+ dt.*v .+ (half_dt_sq).*a

            #  # Prediction
            #  @inbounds @fastmath v .= v .+ half_dt.*a
            #  a .= 0.

            #  # Internal forces -K*d[t+1] stored in global array 'a'
            #  # This is different from matlab code; will change if Nel_ETA is not zero
            #  Threads.@threads for eo in 1:P[1].Nel
                #  ig = P[5][:,:,eo]
                #  a_elem[:,:,eo] = element_computation!(P[2], ig, d[ig], P[7], P[8], P[6][:,:,eo])
            #  end
            #  a = -Conn*a_elem[:]
            #  a[P[4].FltIglobBC] .= 0.

            #  # Absorbing boundaries
            #  @inbounds @fastmath a[P[4].iBcL] .= a[P[4].iBcL] .- P[3].BcLC.*v[P[4].iBcL]
            #  @inbounds @fastmath a[P[4].iBcT] .= a[P[4].iBcT] .- P[3].BcTC.*v[P[4].iBcT]

            #  ###### Fault Boundary Condition: Rate and State #############
            #  @inbounds @fastmath FltVfree .= 2*v[P[4].iFlt] .+ 2*half_dt*a[P[4].iFlt]./P[3].M[P[4].iFlt]
            #  @inbounds @fastmath Vf .= 2*vPre[P[4].iFlt] .+ P[2].Vpl


            #  # Sliprate and NR search
            #  Threads.@threads for j = NFBC:P[1].FltNglob
                #  psi1[j], Vf1[j], tau1[j], psi2[j], Vf2[j], tau2[j] = FBC!(P[1].IDstate, P[3].xLf[j], P[3].Vo[j], P[3].fo[j], P[3].cca[j], P[3].ccb[j], P[3].Seff[j], P[3].tauo[j], P[3].FltZ[j], psi1[j], Vf1[j], tau1[j], psi2[j], Vf2[j], tau2[j], psi[j], Vf[j], FltVfree[j], dt)
            #  end

            #  tau .= tau2 .- P[3].tauo
            #  tau[iFBC] .= 0.
            #  psi .= psi2
            #  @inbounds @fastmath a[P[4].iFlt] .= a[P[4].iFlt] .- P[3].FltB.*tau
            #  ########## End of fault boundary condition ############## 


            #  # Solve for a_new
            #  @inbounds @fastmath a .= a./P[3].M
            
            #  # Correction
            #  @inbounds @fastmath v .= v .+ half_dt*a

            #  v[P[4].FltIglobBC] .= 0.
            #  a[P[4].FltIglobBC] .= 0.

            #  #### Line 861: Omitting P_Ma
            
        #  end # of isolver if loop
        
        #  Vfmax = 2*maximum(v[P[4].iFlt]) .+ P[2].Vpl


        #  #----
        #  # Output variables at different depths for every timestep
        #  # Omitted the part of code from line 871 - 890, because I 
        #  # want to output only certain variables each timestep
        #  # Doing it in separate script
        #  #----


        #  #-----
        #  # Output the variables before and after events
        #  #-----
        #  if Vfmax > 1.01*P[2].Vthres && slipstart == 0
            #  it_s = it_s + 1
            #  delfref = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            #  slipstart = 1
            #  output.tStart[it_s] = output.time_[it]
            #  output.taubefore[:,it_s] = (tau +P[3].tauo)./1e6
            #  vhypo, indx = findmax(2*v[P[4].iFlt] .+ P[2].Vpl)
            #  output.hypo[it_s] = P[3].FltX[indx]
        #  end
        #  if Vfmax < 0.99*P[2].Vthres && slipstart == 1
            #  it_e = it_e + 1
            #  output.delfafter[:,it_e] = 2*d[P[4].iFlt] .+ P[2].Vpl*t .- delfref 
            #  output.tauafter[:,it_e] = (tau + P[3].tauo)./1e6
            #  output.tEnd[it_e] = output.time_[it]
            #  slipstart = 0
        #  end 
        #  #-----
        #  # Output the variables certain timesteps: 2yr interseismic, 1 sec dynamic
        #  #-----
        #  if output.time_[it] > tvsx
            #  ntvsx = ntvsx + 1
            #  idd += 1
            #  output.is_slip[:,ntvsx] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
            #  output.is_slipvel[:,ntvsx] = 2*v[P[4].iFlt] .+ P[2].Vpl
            #  output.is_stress[:,ntvsx] = (tau + P[3].tauo)./1e6
            #  output.index_eq[idd] = 1

            #  tvsx = tvsx + tvsxinc
        #  end

        #  if Vfmax > P[2].Vevne
            #  if idelevne == 0
                #  nevne = nevne + 1
                #  idd += 1
                #  idelevne = 1
                #  tevneb = output.time_[it]
                #  tevne = tevneinc

                #  output.seismic_slip[:,nevne] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
                #  output.seismic_slipvel[:,nevne] = 2*v[P[4].iFlt] .+ P[2].Vpl
                #  output.seismic_stress[:,nevne] = (tau + P[3].tauo)./1e6
                #  output.index_eq[idd] = 2
            #  end

            #  if idelevne == 1 && (output.time_[it] - tevneb) > tevne
                #  nevne = nevne + 1
                #  idd += 1
                
                #  output.seismic_slip[:,nevne] = 2*d[P[4].iFlt] .+ P[2].Vpl*t
                #  output.seismic_slipvel[:,nevne] = 2*v[P[4].iFlt] .+ P[2].Vpl
                #  output.seismic_stress[:,nevne] = (tau + P[3].tauo)./1e6
                #  output.index_eq[idd] = 2
                #  tevne = tevne + tevneinc
            #  end
        #  end

        #  # Output timestep info on screen
        #  if mod(it,500) == 0
            #  @printf("Time (yr) = %1.5g\n\n", t/P[1].yr2sec)
        #  end
        
        #  # Determine quasi-static or dynamic regime based on max-slip velocity
        #  if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
            #  isolver = 1
        #  else
            #  isolver = 2
        #  end
        
        #  output.Vfmax[it] = Vfmax 
       
        #  current_sliprate = 2*v[P[4].iFlt] .+ P[2].Vpl

        #  # Compute next timestep dt
        #  dt = dtevol!(dt , dtmin, P[3].XiLf, P[1].FltNglob, NFBC, current_sliprate, isolver)

    #  end # end of time loop
    
    #  # Remove zeros from preallocated vectors
    #  output.seismic_stress   = output.seismic_stress[:,1:nevne]
    #  output.seismic_slipvel  = output.seismic_slipvel[:,1:nevne]
    #  output.seismic_slip     = output.seismic_slip[:,1:nevne]
    #  output.index_eq         = output.index_eq[1:idd]
    #  output.is_stress        = output.is_stress[:,1:ntvsx]
    #  output.is_slipvel       = output.is_slipvel[:,1:ntvsx]
    #  output.is_slip          = output.is_slip[:,1:ntvsx]
    #  output.tStart           = output.tStart[1:it_s]
    #  output.tEnd             = output.tEnd[1:it_e]
    #  output.taubefore        = output.taubefore[:,1:it_s]
    #  output.tauafter         = output.tauafter[:,1:it_e]
    #  output.delfafter        = output.delfafter[:,1:it_e]
    #  output.hypo             = output.hypo[1:it_s]
    #  output.time_            = output.time_[1:it]
    #  output.Vfmax            = output.Vfmax[1:it]


#  mutable struct results
    #  seismic_stress::Array{Float64,2}
    #  seismic_slipvel::Array{Float64,2}
    #  seismic_slip::Array{Float64,2}
    #  index_eq::Array{Float64}
    #  is_stress::Array{Float64,2}
    #  is_slipvel::Array{Float64,2}
    #  is_slip::Array{Float64,2}
    #  tStart::Array{Float64}
    #  tEnd::Array{Float64}
    #  taubefore::Array{Float64,2}
    #  tauafter::Array{Float64,2}
    #  delfafter::Array{Float64,2}
    #  hypo::Array{Float64}
    #  time_::Array{Float64}
    #  Vfmax::Array{Float64}
#  end

