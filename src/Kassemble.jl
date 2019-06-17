# Assembly of global stiffness matrix as a sparse matrix


function stiffness_assembly(NGLL, NelX, NelY, nglob, dxe, dye, ThickX, ThickY, 
                            rho1, vs1, rho2, vs2, iglob) 
    xgll, wgll, H::SMatrix{NGLL,NGLL,Float64} = GetGLL(NGLL)
    wgll2::SMatrix{NGLL,NGLL,Float64} = wgll*wgll'

    K = SparseMatrixCOO()
    
    ig::Matrix{Int64} = zeros(NGLL,NGLL)  # iterator

    W = material_properties(NelX, NelY,NGLL,dxe, dye, ThickX, ThickY, wgll2, rho1, rho2, vs1, vs2)
    Ke = K_element(W, dxe, dye, NGLL, H, NelX*NelY)
    #  Ksparse = assembley(Ke, iglob, NelX*NelY, nglob)

    # using FEMSparse
    for eo in 1:NelX*NelY
        FEMSparse.assemble_local_matrix!(K, vec(iglob[:,:,eo]),
                                         vec(iglob[:,:,eo]), Ke[:,:,eo])
    end

    return SparseMatrixCSC(K)
end

#  function assembley(Ke, iglob, Nel, nglob)
    #  Ksparse::SparseMatrixCSC{Float64} = spzeros(nglob,nglob) 
    #  for eo in 1:Nel
        #  ig = iglob[:,:,eo]
        #  Ksparse[vec(ig),vec(ig)] += Ke[:,:,eo]
    #  end

    #  Ksparse
#  end


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

function K_element(W, dxe, dye, NGLL, H, Nel)
    # Jacobians
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta

    ww::Matrix{Float64} = zeros(NGLL, NGLL)
    Ke2::Array{Float64,4} = zeros(NGLL,NGLL,NGLL,NGLL)
    Ke::Array{Float64,3} = zeros(NGLL*NGLL,NGLL*NGLL, Nel)
    ig::Matrix{Int64} = zeros(NGLL,NGLL)  # iterator
    
    term1::Float64 = 0.; term2::Float64 = 0.
    del = Matrix{Float64}(I,NGLL,NGLL)  # identity matrix
    
    #  @threads for tid in 1:nthreads()
        #  len = div(NelX*NelY,nthreads())
        #  domain = ((tid-1)*len +1):tid*len

        @inbounds for eo in 1:Nel
            Ke2 .= 0.

            ww = W[:,:,eo]
            term1 = 0.; term2 = 0.
            @inbounds for i in 1:NGLL, j in 1:NGLL
                term1 = 0; term2 = 0
                @inbounds for k in 1:NGLL, l in 1:NGLL
                    term1 = 0; term2 = 0
                    @simd for p in 1:NGLL
                        term1 += del[i,k]*ww[k,p]*(jac/dy_deta^2)*H[j,p]*H[l,p]
                        term2 += del[j,l]*ww[p,j]*(jac/dx_dxi^2)*H[i,p]*H[k,p]
                    end
                    @inbounds Ke2[i,j,k,l] = term1 + term2
                end
            end
            Ke[:,:,eo] = reshape(Ke2,NGLL*NGLL,NGLL*NGLL)
        end
    #  end

    Ke
end