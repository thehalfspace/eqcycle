####################################
#   NEWTON RHAPSON SEARCH METHOD
####################################

# Fault Boundary function
function FBC!(IDstate::Int, P::params_farray, NFBC::Int, FltNglob::Int, psi1::Array{Float64}, Vf1::Array{Float64}, tau1::Array{Float64}, psi2::Array{Float64}, Vf2::Array{Float64}, tau2::Array{Float64}, psi::Array{Float64}, Vf::Array{Float64}, FltVfree::Array{Float64}, dt::Float64)

    tauNR::Vector{Float64} = zeros(FltNglob)

   for j = NFBC:FltNglob 

        psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-5, IDstate)

        Vf1[j], tau1[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],
                                    tauNR[j], P.tauo[j], psi1[j], P.FltZ[j], FltVfree[j])
    
        if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
            
            println("Fault Location = ", j)
            println(" Vf = ", Vf[j])
            println(" tau1 = ", tau1[j])

            println("psi =", psi[j])
            println("psi1 =", psi1[j])
            # Save simulation results
            #filename = string(dir, "/data", name, "nrfail.jld2")
            #@save filename results(Stress,SlipVel, Slip, time_) 
            @error("NR SEARCH FAILED!")
            return
        end
        
        psi2[j] = IDS2!(P.xLf[j], P.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], IDstate)
        
        # NRsearch 2nd loop
        Vf2[j], tau2[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],
                                  tau1[j], P.tauo[j], psi2[j], P.FltZ[j], FltVfree[j])

    end

    return psi1, Vf1, tau1, psi2, Vf2, tau2
end


# Newton Rhapson search method
function NRsearch!(fo::Float64, Vo::Float64, cca::Float64, ccb::Float64, Seff::Float64, tau::Float64, tauo::Float64, psi::Float64, FltZ::Float64, FltVfree::Float64)

    Vw = 1e10
    fact = 1 + (Vo/Vw)*exp(-psi)

    # NR search point by point for tau if Vf < Vlimit
    eps = 0.001*cca*Seff
    k = 0
    delta = Inf

    while abs(delta) > eps
        fa = fact*tau/(Seff*cca)
        help = -(fo + ccb*psi)/cca

        help1 = exp(help + fa)
        help2 = exp(help - fa)

        Vf = Vo*(help1 - help2)

        Vfprime = fact*(Vo/(cca*Seff))*(help1 + help2)

        delta = (FltZ*FltVfree - FltZ*Vf + tauo - tau)/(1 + FltZ*Vfprime)

        tau = tau + delta
        k = k + 1

        if abs(delta) > 1e10 || k == 1000
            println("k = ", k)
            # Save simulation results
            #filename = string(dir, "/data", name, "nrfail.jld2")
            #@save filename 
            @error("NR search fails to converge")
        end
    end

    fa = fact*tau/(Seff*cca)
    
    help = -(fo + ccb*psi)/cca

    help1 = exp(help + fa)
    help2 = exp(help - fa)

    Vf = Vo*(help1 - help2)

    return Vf, tau
end
