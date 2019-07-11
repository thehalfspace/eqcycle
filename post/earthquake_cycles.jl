#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################

using StatsBase
using PyPlot

PyPlot.matplotlib.rc("patch.force_edgecolor=true")

# get index of start of rupture
function get_index(seismic_stress, taubefore)

    index_start = zeros(Int, length(taubefore[1,:]))
    for i in 2:length(taubefore[1,:])
        index_start[i] = findall(O.seismic_stress[1921,:] .== O.taubefore[1921,i])[1]
    end

    index_start
end

#.................................................
# Compute the final Coseismic slip for each event
#.................................................
function Coslip(S, Slip, SlipVel, Stress, time_=zeros(1000000))
    Vfmax = maximum(SlipVel, dims = 1)[:]

    delfafter::Array{Float64,2} = zeros(size(Slip))
    tStart::Array{Float64} = zeros(size(Slip[1,:]))
    tEnd::Array{Float64} = zeros(size(Slip[1,:]))

    taubefore::Array{Float64,2} = zeros(size(Slip))
    tauafter::Array{Float64,2} = zeros(size(Slip))
    
    hypo::Array{Float64} =  zeros(size(Slip[1,:]))   # Hypocenter
    vhypo::Array{Float64} = zeros(size(Slip[1,:]))   # Velocity at hypocenter

    Vthres = 0.001 # event threshold
    slipstart = 0
    it = 1; it2 = 1
    delfref = zeros(size(Slip[:,1]))

    for i = 1:length(Slip[1,:])

        # Start of each event
        if Vfmax[i] > 1.01*Vthres && slipstart == 0
            delfref = Slip[:,i]
            slipstart = 1
            tStart[it2] = time_[i]
            
            taubefore[:,it2] = Stress[:,i]
            vhypo[it2], indx = findmax(SlipVel[:,i])

            hypo[it2] = S.FltX[indx]

            it2 = it2+1
        end

        # End of each event
        if Vfmax[i] < 0.99*Vthres && slipstart == 1
            delfafter[:,it] = Slip[:,i] - delfref
            tauafter[:,it] = Stress[:,i]
            tEnd[it] = time_[i]
            slipstart = 0
            it = it + 1
        end
    end

    return delfafter[:,1:it-1], (taubefore-tauafter)[:,1:it-1], tStart[1:it2-1], tEnd[1:it-1], vhypo[1:it2-1], hypo[1:it2-1]
end

#..........................................................
# Compute the moment magnitude:
#       Assumed the rupture area to be square; the rupture
#       dimension along depth is the same as the rupture
#       dimension perpendicular to the plane
#..........................................................
function moment_magnitude_new(mu, P3, FltX, delfafter, stressdrops ,time_)
    # Final coseismic slip of each earthquake
    #  delfafter, stressdrops = Coslip(S, Slip, SlipVel, Stress, time_)

    iter = length(delfafter[1,:])
    seismic_moment = zeros(iter)
    fault_slip = zeros(iter)
    temp_sigma = 0
    iter2 = 1 

    del_sigma = zeros(iter)
    
    dx = diff(FltX)

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0; temp_slip = 0

        for j = 1:P3.FltNglob
            if delfafter[j,i] >= slip_thres
                area = area + delfafter[j,i]*dx[j-1]
                zdim = zdim + dx[j-1]
                temp_slip = temp_slip + delfafter[j,i]

                # Avg. stress drops along rupture area
                temp_sigma = temp_sigma + stressdrops[j,i]*dx[j-1]
            end
        end
        
        seismic_moment[i] = mu*area*zdim
        del_sigma[i] = temp_sigma/zdim
        fault_slip[i] = temp_slip/zdim


    end
    #  seismic_moment = filter!(x->x!=0, seismic_moment)
    #  del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7

    return Mw, del_sigma, fault_slip
end


#--------------------------
# Earthquake scaling laws
#--------------------------
function scaling(a,b)

    fig = PyPlot.figure(figsize=(12,9))
    ax = fig.add_subplot(111)

    ax.plot(a,b, "k.", markersize=20)
    ax.set_xlabel("Slip (m)")
    ax.set_ylabel("Rupture Length (m)")
    #  ax.set_yscale("log")
    #  ax.set_xscale("log")
    ax.set_title("Rupture length vs Slip")
    #  ax.legend(loc="upper right")
    show()

    figname = string(path, "scaling1.png")
    fig.savefig(figname, dpi = 300)
end

