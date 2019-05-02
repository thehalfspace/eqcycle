################################################
#                                              
#   SOLVE FOR DISPLACEMENT USING PRECONDITIONED 
#           CONJUGATE GRADIENT METHOD          
#                                              
################################################

using IterativeSolvers
#  using AlgebraicMultigrid

# trying matrix free operatiors
struct MatrixOperator4
    P::params_float
    iglob::Array{Int64, 3}
    H::Array{Float64,2}
    Ht::Array{Float64,2}
    W::Array{Float64,3}
    Nel::Int64
    FltNI::Array{Int64}
    size::Tuple{Int, Int}
end


Base.eltype(A::MatrixOperator4) = Float64
Base.size(A::MatrixOperator4) = A.size
Base.size(A::MatrixOperator4,i::Int) = A.size[i]

function conj_grad(P, F)
    
    W = P[6]
    P2 = P[2];
    iglob = P[5];
    H = P[7]
    Ht = P[8]
    Nel = P[1].Nel
    nglob = P[1].nglob
    FltNI = P[4].FltNI
    p = P[3].diagKnew

    #  a_local::Array{Float64} = zeros(nglob)
    
    A = MatrixOperator4(P2,iglob, H, Ht, W, Nel, FltNI, (nglob, nglob))
    #  p = aspreconditioner(A)
    #  LinearAlgebra.mul!(a_local,A,dnew)
    #  Fnew = -a_local[P[4].FltNI]  # b 
   
    #  return A,Fnew
    U = cg(A,F)
    return  U[FltNI]


end

function LinearAlgebra.mul!(C, A::MatrixOperator4, b)
    #  Threads.@threads for tid in 1:Threads.nthreads()
        #  len = div(A.Nel, Threads.nthreads())
        #  domain = ((tid-1)*len + 1):tid*len

        @inbounds @simd for eo in 1:A.Nel
        ig = A.iglob[:,:,eo]
        Wlocal = A.W[:,:,eo]
        locall = b[ig]
        C[ig] = C[ig] +  A.P.coefint1*A.H*(Wlocal.*(A.Ht*locall)) + 
                    A.P.coefint2*(Wlocal.*(locall*A.H))*A.Ht
        end
    #  end
    return C
end
LinearAlgebra.:*(A::MatrixOperator4,B::AbstractVector) = (C = similar(B); LinearAlgebra.mul!(C,A,B))

