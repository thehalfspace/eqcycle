############
# Try assembling stiffness
###########
function assemble(NGLL, iglob, jac)

    Nsize = 2
    LX::Int = 24e3*Nsize  # depth dimension of rectangular domain
    LY::Int = 15e3*Nsize # off fault dimenstion of rectangular domain

    res = 1
    NelX::Int = 15*Nsize*res # no. of elements in x
    NelY::Int = 10*Nsize*res # no. of elements in y
    Nel::Int = NelX*NelY # Total no. of elements

    dxe::Float64 = LX/NelX #	Size of one element along X
    dye::Float64 = LY/NelY #	Size of one element along Y

    ThickX = 48e3
    ThickY = 0
    # default
    rho1::Float64 = 2670
    vs1::Float64 = 3464

    rho2::Float64 = 2500
    vs2::Float64 = 0.6*vs1

    xgll, wgll, H = GetGLL(NGLL)
    wgll2 = wgll*wgll';
    W::Array{Float64,3} = zeros(NGLL, NGLL, Nel)

    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    
    for ey = 1:NelY
        for ex = 1:NelX

            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]

            # Properties of heterogeneous medium
            if ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)
                rho[:,:] .= rho2
                mu[:,:] .= rho2*vs2^2
            else
                rho[:,:] .= rho1
                mu[:,:] .= rho1*vs1^2
            end

            # Local contributions to the stiffness matrix
            KK[:,:,eo] .= wgll2.*mu*jac;

        end
    end

    return KK
end
