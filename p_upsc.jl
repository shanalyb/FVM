function upscaling(k::Array{Float64,2}, N::Int64)
    k_new = zeros(2*N-1, 1)
    for i = 1:N
        k_new[2*i-1] = k[i]
    end
    for i = 1:N-1
        k_new[2*i] = abs((k[i+1] + k[i]) / 2)
    end
    k_new
end
