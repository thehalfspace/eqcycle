#####################################
# DAMAGE EVOLUTION IN TIME
#####################################

function damageEvol!(damage_idx, W)
    W[:,:,damage_idx] .= W[:,:,damage_idx]*0.95
end
