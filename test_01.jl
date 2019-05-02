using LinearAlgebra
using LinearMaps
using IterativeSolvers

n = 100
#  A = rand(n, n)
B = rand(n)

struct someMatrix
    size::Tuple{Int,Int}
end

Base.eltype(A::someMatrix) = Float64
Base.size(A::someMatrix) = A.size
Base.size(A::someMatrix,i::Int) = A.size[i]

A = someMatrix((length(B),length(B)))

function LinearAlgebra.mul!(C,A::someMatrix,B)
    for i in 2:length(B)-1
        C[i] = B[i-1] - 2B[i] + B[i+1]
    end
    C[1] = -2B[1] + B[2]
    C[end] = B[end-1] - 2B[end]
    C
end
Base.:*(A::someMatrix,B::AbstractVector) = (C = similar(B); mul!(C,A,B))


U = gmres(A,B)
norm(A*U - B)
