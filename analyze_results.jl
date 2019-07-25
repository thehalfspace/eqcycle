############################################
### ANALYZE RESULTS FROM THE OUTPUT FILES
############################################

include("output.jl")

include("post/earthquake_cycles.jl")
include("post/plots.jl")
include("post/cumulative_slip.jl")

# path to save files
global path = "$(@__DIR__)/plots/test01/"

# Deserialize the output
using Serialization
open("data/dc8_b20.out") do f
    global O, sim_time, P, S
    O = deserialize(f)
    sim_time = deserialize(f)
    P = deserialize(f)
end


P1 = P[1]
P2 = P[2]
P3 = P[3]
P4 = P[4]

FltX =P3.FltX

delfsec = O.seismic_slip
delf5yr = O.is_slip
delfafter = O.delfafter
stressdrops = O.taubefore-O.tauafter
seismic_stress = O.seismic_stress
seismic_slipvel = O.seismic_slipvel
seismic_slip = O.seismic_slip
index_eq = O.index_eq
is_stress = O.is_stress
is_slipvel = O.is_slipvel
is_slip = O.is_slip
tStart = O.tStart
tEnd = O.tEnd
delfafter = O.delfafter
hypo = O.hypo
time_ = O.time_
Vfmax = O.Vfmax
yr2sec = P1.yr2sec

# time-index of start of rupture
start_index = get_index(seismic_stress, O.taubefore)

rho1 = 2670
vs1 = 3464
rho2 = 2500
vs2 = 0.6*vs1
mu = rho2*vs2^2

Mw, del_sigma, fault_slip = moment_magnitude_new(mu, P1, FltX, delfafter, stressdrops, O.time_);


#-----------------------------
# test new code snippets here
#-----------------------------
using Plots
# animate rupture tip
function anim_rupture(stress, FltX)
    
    anim = @animate for i = 1:length(stress[1,:])
        Plots.plot(stress[:,i], yaxis=("Depth (km)", -FltX./1e3, :flip), leg=false)

    end

    gif(anim, "$(@__DIR__)/plots/temp.gif", fps=5)

end

