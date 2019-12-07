
# Material properties for a narrow rectangular damaged zone of half-thickness ThickY and depth ThickX 
function material_properties(NelX, NelY,NGLL, dxe, dye, ThickX, ThickY, wgll2, rho1, rho2, vs1, vs2)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            
            # Properties of heterogeneous medium
            if ex*dxe >= ThickX && (dye <= ey*dye <= ThickY)
                mu .= rho2*vs2^2
            else
                mu .= rho1*vs1^2
            end
            W[:,:,eo] = wgll2.*mu
        end
    end
    W
end



# Material properties for a trapezium damage zone
# Linear function: shape of trapezoid
function line(x,y)
    P1 = [0 1.5e3]
    P2 = [-8e3 150]

    f = (y - P2[2]) - ((P1[2]-P2[2])/(P1[1]-P2[1]))*(x - P2[1])

    return f
end

# Set up trapezoidal rigidity
function rigid(x,y)
    # Rigidity: host rock and fault zone
    rho1::Float64 = 2670
    vs1::Float64 = 3464
    
    rho2 = 0.60*rho1
    vs2 = 0.60*vs1
    rho3 = 0.8*rho1
    vs3 = 0.8*vs1
    
    rhoglob::Array{Float64} = zeros(length(x))
    vsglob::Array{Float64} = zeros(length(x))

    for i = 1:length(x)
        if x[i] > -8e3
            if line(x[i],y[i]) < 0
                rhoglob[i] = rho3
                vsglob[i] = vs3
            else
                rhoglob[i] = rho1
                vsglob[i] = vs1
            end
        else
            rhoglob[i] = rho1
            vsglob[i] = vs1
        end

    end

    for i = 1:length(x)
        if y[i]< 150
            rhoglob[i] = rho2
            vsglob[i] = vs2
        end
    end

    
    return rhoglob, vsglob
end

function mat_trap(NelX, NelY, NGLL, iglob, M, dxe, dye, x,y, wgll2)
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    rhoglob, vsglob = rigid(x,y)
    muglob = rhoglob.*(vsglob.^2)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]
            
            mu[:,:] = muglob[ig]
            rho[:,:] = rhoglob[ig]
            
            W[:,:,eo] = wgll2.*mu
            M[ig] .+= wgll2.*rho*jac
        end
    end
    return M,W
end

# Setup for Gaussian rigidity
function mat_gauss(NelX, NelY, NGLL, iglob, M, dxe, dye, x,y, wgll2)
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    
    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    rhoglob, vsglob = rigid_gauss(x,y)
    muglob = rhoglob.*(vsglob.^2)
    
    @inbounds for ey in 1:NelY
        @inbounds for ex in 1:NelX
            eo = (ey-1)*NelX + ex
            ig = iglob[:,:,eo]
            
            mu[:,:] = muglob[ig]
            rho[:,:] = rhoglob[ig]
            
            W[:,:,eo] = wgll2.*mu
            M[ig] .+= wgll2.*rho*jac
        end
    end
    return M,W
end

# Set up gaussian rigidity
function rigid_gauss(x,y)

    # Rigidity: host rock and fault zone
    rho1::Float64 = 2670
    vs1::Float64 = 3464
    
    rho2 = 0.65*rho1
    vs2 = 0.65*vs1
    muhost = rho1*vs1^2
    mufz = rho2*vs2^2

    LX = 48.0e3
    ThickX = 12.0e3
    ThickY = 1.0e3

    # Gaussian fault zone mean and std
    meanx = 0
    meany = 0
    sigx = (LX - ThickX)/3
    sigy = ThickY/3

    muglob = (mufz-muhost)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ muhost
    
    vsglob = (vs2-vs1)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ vs1
    rhoglob = (rho2-rho1)*exp.(-(gauss(x, meanx, sigx) .+
                        gauss(y, meany, sigy))) .+ rho1
    
    return vsglob, rhoglob
end


# Gaussian function
function gauss(x, mu, sigma)
    return ((x .- mu)./(2*sigma)).^2
end
