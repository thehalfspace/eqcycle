n = 1000000


arr = zeros(n)

function single(n,arr)
    for i in 1:n
        arr[i] = i^2
    end
end

function multi(n,arr)
    @threads for tid in 1:nthreads()
        len = div(n, nthreads())
        domain = ((tid-1)*len +1):tid*len
        @inbounds for i in domain
            arr[i] = i^2
        end
    end
end
