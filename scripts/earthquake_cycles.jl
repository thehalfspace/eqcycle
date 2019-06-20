#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################

using StatsBase
using PyPlot

PyPlot.matplotlib[:rc]("patch.force_edgecolor=true")

#-------------------------------
# Compute hypocenter locations
#------------------------------
function plotHypo(hypo)  #S, Slip, SlipVel, Stress, time_)

    # Plot hypocenter
    hist = fit(Histogram, -hypo./1e3, closed=:right, nbins=10)

    fig = PyPlot.figure(figsize=(5,7))
    ax = fig[:add_subplot](111)

    ax[:barh](hist.edges[1][1:end-1], hist.weights)
    #  ax[:plot](collect(1:80), -8*ones(80), "k--", label="Fault Zone Depth")
    ax[:set_xlabel]("Number of Earthquakes")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Hypocenter Location")
    ax[:set_ylim]([0, 24])
    ax.invert_yaxis()
    fig.tight_layout()
    show()

    figname = string(path, "hypo.png")
    fig[:savefig](figname, dpi = 300)

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
function moment_magnitude(P, S, delfafter, stressdrops ,time_)
    # Final coseismic slip of each earthquake
    #  delfafter, stressdrops = Coslip(S, Slip, SlipVel, Stress, time_)

    iter = length(delfafter[1,:])
    mu = P.rho1*P.vs1^2
    seismic_moment = zeros(iter)
    rupture_length = zeros(iter)
    fault_slip = zeros(iter)
    temp_sigma = 0
    iter2 = 1 

    del_sigma = zeros(iter)
    
    dx = diff(S.FltX)

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0; temp_slip=0;

        for j = 1:P.FltNglob
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
        rupture_length[i] = zdim

    end
    #  seismic_moment = filter!(x->x!=0, seismic_moment)
    #  del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7

    return Mw, del_sigma, fault_slip, rupture_length
end

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
    ax = fig[:add_subplot](111)

    ax[:plot](a,b, "k.", markersize=20)
    ax[:set_xlabel]("Slip (m)")
    ax[:set_ylabel]("Rupture Length (m)")
    #  ax[:set_yscale]("log")
    #  ax[:set_xscale]("log")
    ax[:set_title]("Rupture length vs Slip")
    #  ax[:legend](loc="upper right")
    show()

    figname = string(path, "scaling1.png")
    fig[:savefig](figname, dpi = 300)
end


#...........
# Plot MFD
#...........
function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 9)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(8,7))
    ax = fig[:add_subplot](111)

    #  ax[:plot](hist.edges[1][1:end-1], hist.weights, ".", label="Non-cumulative")
    ax[:plot](hist.edges[1][1:end-1], cum, "k.--", markersize=20) #, label="Cumulative")
    #  ax[:set_xlabel]("Moment Magnitude (Mw)")
    #  ax[:set_ylabel]("Number of Earthquakes")
    ax[:set_yscale]("log")
    #  ax[:set_title]("Magnitude-frequency distribution")
    #  ax[:set_xlim]([2, 7])
    ax[:set_ylim]([1, 30])
    #  ax[:legend](loc="upper right")
    show()

    figname = string(path, "mfd.png")
    fig[:savefig](figname, dpi = 300)
end

function MwPlot2(Mw1, Mw2, Mw3, Mw4, Mw5, Mw6)

    hist1 = fit(Histogram, Mw1, nbins = 6)
    hist2 = fit(Histogram, Mw2, nbins = 6)
    hist3 = fit(Histogram, Mw3, nbins = 6)
    hist4 = fit(Histogram, Mw4, nbins = 6)
    hist5 = fit(Histogram, Mw5, nbins = 6)
    hist6 = fit(Histogram, Mw5, nbins = 6)

    # Cumulative
    cum1 = cumsum(hist1.weights[end:-1:1])[end:-1:1]
    cum2 = cumsum(hist2.weights[end:-1:1])[end:-1:1]
    cum3 = cumsum(hist3.weights[end:-1:1])[end:-1:1]
    cum4 = cumsum(hist4.weights[end:-1:1])[end:-1:1]
    cum5 = cumsum(hist5.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(12,9))
    ax = fig[:add_subplot](111)

    ax[:plot](hist1.edges[1][1:end-1], cum1, ".-", markersize=10, label="400 m")
    ax[:plot](hist2.edges[1][1:end-1], cum2, ".-", markersize=10, label="320 m")
    ax[:plot](hist3.edges[1][1:end-1], cum3, ".-", markersize=10, label="250 m")
    ax[:plot](hist4.edges[1][1:end-1], cum4, ".-", markersize=10, label="160 m")
    ax[:plot](hist5.edges[1][1:end-1], cum5, ".-", markersize=10, label="100 m")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel]("Number of Earthquakes")
    ax[:set_yscale]("log")
    ax[:set_title]("Magnitude-frequency distribution")
    #  ax[:set_xlim]([1.7, 7])
    ax[:set_ylim]([1e-1, 1e3])
    ax[:legend](loc="upper right")
    show()

    figname = string(path, "mfdcomp.png")
    fig[:savefig](figname, dpi = 300)
end

#.................................
# Plot earthquake catalog
# (Earthquake magnitude with time)
#.................................
function eq_catalog(Mw, t_catalog, yr2sec)

    fig = PyPlot.figure(figsize=(12,9))
    ax = fig[:add_subplot](111)

    ax[:scatter](t_catalog./yr2sec, Mw, s= 30, marker=".")
    ax[:set_xlabel]("Time (yrs)")
    ax[:set_ylabel]("Moment Magnitude (Mw)")
    ax[:set_title]("Earthquake Catalogue")
    show()

    figname = string(path, "catalogue.pdf")
    fig[:savefig](figname, dpi = 300)
end
