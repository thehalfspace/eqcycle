########################################
#   ISOLATE EACH EARTHQUAKE STUFF
########################################

function event_indx2(O, P)
    start_indx = zeros(size(O.tStart))
    end_indx = zeros(size(O.tStart))
    evno = zeros(size(O.tStart))

    # indx = 2 for aseismic, 1 for seismic
    indx::Array{Float64} = zeros(length(O.time_))
    
    yr2sec = 365*24*60*60 
    tvsx = 2*P.yr2sec
    tvsxinc = tvsx

    tevneinc = 1
    Vevne = 0.001

    ntvsx = 0
    n3 = 0
    nevne = 0
    idelevne = 0
    tevneb = 0
    tevne = 0
    
    it = length(O.time_)
  
    for i = 1:it
        if O.time_[i] > tvsx
            ntvsx = ntvsx + 1
            n3 = n3 + 1
            indx[n3] = 2     # interseismic
            tvsx = tvsx +tvsxinc
        end
        
        if O.Vfmax[i] > Vevne 
            if idelevne == 0
                nevne = nevne + 1
                idelevne = 1
                tevneb = O.time_[i]
                tevne = tevneinc
                n3 = n3 + 1
                indx[n3] = 1     # seismic
            end

            if idelevne == 1 && (O.time_[i] - tevneb) > tevne
                nevne = nevne + 1
                n3 = n3 + 1
                indx[n3] = 1
                #Vfsec[:,nevne] = 2*v[S.iFlt] .+ P.Vpl
                #Tausec[:,nevne] = (tau + tauo)./1e6
                tevne = tevne + tevneinc
            end

        else
            idelevne = 0
        end
    end 

    return indx[1:n3]
end


# get stress/sliprate for one event
function count_event(O, indx)
    counter = zeros(size(O.tStart))
    j = 1; iterator = 1
    for i in indx
        if i == 2
            counter[j] += 1
            iterator = 1
        else
            if iterator == 1
                println(i)
                j += 1
                iterator = 0
            end
        end
    end
    return counter
end


function slipvel_event(O, indx)

    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(delf5yr, FltX/1e3, color="royalblue", lw=1, alpha=1.0)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Cumulative Slip History")
    ax.set_ylim([-24, 0])
    ax.set_xlim([1,15])  #[0,maximum(delf5yr)])
    show()
    
    figname = string(path, "slipvel_event.pdf")
    #  fig.savefig(figname, dpi = 300)

end
