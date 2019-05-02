
# IDstate functions
function IDS!(xLf::Float64, Vo::Float64, psi::Float64, dt::Float64, Vf::Float64, cnd::Float64, IDstate = 2)
    #= compute slip-rates on fault based on different
       formulations =#

    if IDstate == 1
        psi1 = psi + dt*((Vo./xLf).*exp(-psi) - abs(Vf)./xLf)

    elseif IDstate == 2
        VdtL = abs(Vf)*dt/xLf
        if VdtL < cnd
            psi1 = log( exp(psi-VdtL) + Vo*dt/xLf -
                        0.5*Vo*abs(Vf)*dt*dt/(xLf^2))
        else
            psi1 = log(exp(psi-VdtL) + (Vo/abs(Vf))*(1-exp(-VdtL)))
        end

    elseif IDstate == 3
        psi1 = exp(-abs(Vf)*dt/xLf) * log(abs(Vf)/Vo) + 
        exp(-abs(Vf)*dt/xLf)*psi + log(Vo/abs(Vf))

        if ~any(imag(psi1)) == 0
            return
        end
    end

    return psi1

end

# On fault slip rates
function IDS2!(xLf::Float64, Vo::Float64, psi::Float64, psi1::Float64, dt::Float64, Vf::Float64, Vf1::Float64, IDstate = 2)
            
    if IDstate == 1
        psi2 = psi + 0.5*dt*( (Vo/xLf)*exp(-psi) - abs(Vf)/xLf 
                                + (Vo/xLf)*exp(-psi1) - abs(Vf1)/xLf )

    elseif IDstate == 2
        VdtL = 0.5*abs(Vf1 + Vf)*dt/xLf

        if VdtL < 1e-6
            psi2 = log( exp(psi-VdtL) + Vo*dt/xLf -
                            0.5*Vo*0.5*abs(Vf1 + Vf)*dt*dt/(xLf^2))
        else
            psi2 = log(exp(psi-VdtL) + 
                            (Vo/(0.5*abs(Vf + Vf1)))*(1-exp(-VdtL)))
        end

    elseif IDstate == 3
        psi2 = exp(-0.5*abs(Vf + Vf1)*dt/xLf) * log(0.5*abs(Vf + Vf1)/Vo) + 
                        exp(-0.5*abs(Vf + Vf1)*dt/xLf)*psi 
                        + log(Vo/(-0.5*abs(Vf + Vf1)) )
    end

    return psi2
end

# Slip rates on fault for quasi-static regime
function slrFunc!(P::params_farray, NFBC::Int, FltNglob::Int, psi::Array{Float64}, psi1::Array{Float64}, Vf::Array{Float64}, Vf1::Array{Float64}, IDstate::Int, tau1::Array{Float64}, dt::Float64)
    
    Threads.@threads for tid in 1:Threads.nthreads()
        len = div(length(NFBC:FltNglob), Threads.nthreads())
        domain = ((tid-1)*len + NFBC):(tid*len + NFBC - 1)

        @inbounds @simd for j in domain
            psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)
            Vf1[j] = P.Vo[j]*(exp(-(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j] + 
                                  (tau1[j] + P.tauo[j])/(P.Seff[j]*P.cca[j])) 
                              - exp(-(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j] - 
                                    (tau1[j] + P.tauo[j])/(P.Seff[j]*P.cca[j])))
        end
    end
    j = FltNglob
    psi1[j] = IDS!(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)
    fa = (tau1[j] + P.tauo[j])/(P.Seff[j]*P.cca[j])
    help = -(P.fo[j] + P.ccb[j]*psi1[j])/P.cca[j]
    help1 = exp(help + fa)
    help2 = exp(help - fa)
    Vf1[j] = P.Vo[j]*(help1 - help2) 
    return psi1, Vf1
end
