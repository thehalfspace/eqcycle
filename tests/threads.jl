using LinearAlgebra
using Base.Threads

function run(N)
    a::Array{Float64} = rand(5)
    b::Array{Float64,3} = rand(5,5,N)
    c::Array{Float64} = zeros(N)

    thread1(a,b,c, N)

end

function thread1(a,b,c, N)
    for tid in nthreads()
        len = div(N, nthreads())
        domain = ((tid-1)*len + 1):tid*len

        @inbounds @simd for eo in domain
            
            c[eo] = norm(b[:,:,eo]*a)
        end
    end
end
