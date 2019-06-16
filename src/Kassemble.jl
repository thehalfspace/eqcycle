# Assembly of global stiffness matrix as a sparse matrix
using StaticArrays
function stiffness_assembly(NGLL, NelX, NelY, nglob, dxe, dye, ThickX, ThickY, 
                            rho1, vs1, rho2, vs2, iglob) 
    xgll, wgll, H::SMatrix{NGLL,NGLL,Float64} = GetGLL(NGLL)
    wgll2::SMatrix{NGLL,NGLL,Float64} = wgll*wgll'

    # Jacobians
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    
    #  rho::Matrix{Float64} = zeros(NGLL, NGLL)

    mu::Matrix{Float64} = zeros(NGLL, NGLL)
    W::Array{Float64,3} = zeros(NGLL, NGLL, NelX*NelY)
    ww::Matrix{Float64} = zeros(NGLL, NGLL)

    Ksparse::SparseMatrixCSC{Float64} = spzeros(nglob,nglob) 
    Ke2::Array{Float64,4} = zeros(NGLL,NGLL,NGLL,NGLL)
    Ke::Array{Float64,3} = zeros(NGLL*NGLL,NGLL*NGLL, NelX*NelY)
            
    term1::Float64 = 0.; term2::Float64 = 0.
    del = Matrix{Float64}(I,NGLL,NGLL)  # identity matrix
    
    ig::Matrix{Int64} = zeros(NGLL,NGLL)  # iterator

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

    @threads for tid in 1:nthreads()
        len = div(NelX*NelY,nthreads())
        domain = ((tid-1)*len +1):tid*len

        @inbounds for eo in domain #1:NelX*NelY

            ig = iglob[:,:,eo]

            Ke2 .= 0.

            ww = W[:,:,eo]
            #  W = wgll2.*mu
            term1 = 0.; term2 = 0.
            @inbounds for i in 1:NGLL
                @inbounds for j in 1:NGLL
                    term1 = 0; term2 = 0
                    @inbounds for k in 1:NGLL
                        @inbounds for l in 1:NGLL
                            term1 = 0; term2 = 0
                            @inbounds for p in 1:NGLL
                            term1 += del[i,k]*ww[k,p]*(jac/dy_deta^2)*H[j,p]*H[l,p]
                            term2 += del[j,l]*ww[p,j]*(jac/dx_dxi^2)*H[i,p]*H[k,p]
                            end
                            Ke2[i,j,k,l] = term1 + term2
                        end
                    end
                end
            end

            Ke[:,:,eo] = reshape(Ke2,NGLL*NGLL,NGLL*NGLL)
        end
    end


    @inbounds @fastmath for eo in 1:NelX*NelY
        ig = iglob[:,:,eo]
        Ksparse[ig[:],ig[:]] += Ke[:,:,eo]
    end

    Ksparse
end
