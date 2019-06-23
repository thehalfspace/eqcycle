####################################
#   NEWTON RHAPSON SEARCH METHOD
####################################

# Fault Boundary function
#  sin1(x::Float64) = ccall(:sin, Float64, (Float64,), x)
#  cos1(x::Float64) = ccall(:cos, Float64, (Float64,), x)
exp1(x::Float64) = ccall(:exp, Float64, (Float64,), x)
log1(x::Float64) = ccall(:log, Float64, (Float64,), x)

function FBC!(IDstate, P::params_farray, NFBC, FltNglob, psi1, Vf1, tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)

        @threads for tid in 1:nthreads()
            len = div(length(NFBC:FltNglob), nthreads())
            domain = ((tid-1)*len + NFBC):(tid*len + NFBC - 1)

            @inbounds for j in domain

                psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-5, IDstate)
                Vf1[j], tau1[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j], 0. , P.tauo[j], psi1[j], P.FltZ[j], FltVfree[j])
    
                if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
                    println("Fault Location = ", j)
                    println("Fault Location = ", j)
                    println(" Vf = ", Vf[j])
                    println(" tau1 = ", tau1[j])
                    @error("NR SEARCH FAILED!")
                    return
                end
        
                psi2[j] = IDS2!(P.xLf[j], P.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], IDstate)
                Vf2[j], tau2[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],tau1[j], P.tauo[j], psi2[j], P.FltZ[j], FltVfree[j])
            end
            
        end

        j = FltNglob
        psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-5, IDstate)
        Vf1[j], tau1[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],0., P.tauo[j], psi1[j], P.FltZ[j], FltVfree[j])

        if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
            println("Fault Location = ", j)
            println("Fault Location = ", j)
            println(" Vf = ", Vf[j])
            println(" tau1 = ", tau1[j])
            @error("NR SEARCH FAILED!")
            return
        end

        psi2[j] = IDS2!(P.xLf[j], P.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], IDstate)
        Vf2[j], tau2[j] = NRsearch!(P.fo[j], P.Vo[j], P.cca[j], P.ccb[j], P.Seff[j],tau1[j], P.tauo[j], psi2[j], P.FltZ[j], FltVfree[j])


    return psi1, Vf1, tau1, psi2, Vf2, tau2
end


# Newton Rhapson search method
function NRsearch!(fo::Float64, Vo::Float64, cca::Float64, ccb::Float64, Seff::Float64, tau::Float64, tauo::Float64, psi::Float64, FltZ::Float64, FltVfree::Float64)

    Vw = 1e10
    fact = 1 + (Vo/Vw)*exp1(-psi)

    # NR search point by point for tau if Vf < Vlimit
    eps = 0.001*cca*Seff
    k = 0
    delta = Inf

    while abs(delta) > eps
        fa = fact*tau/(Seff*cca)
        help = -(fo + ccb*psi)/cca

        help1 = exp1(help + fa)
        help2 = exp1(help - fa)

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

    help1 = exp1(help + fa)
    help2 = exp1(help - fa)

    Vf = Vo*(help1 - help2)

    return Vf, tau
end
