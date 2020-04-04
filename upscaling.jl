using StatsBase
function geometric_average(k_old::Array, x::Int64, y::Int64)
    N_x = size(k_old, 1) ÷ x
    N_y = size(k_old, 2) ÷ y
    k = zeros(N_x, N_y)
    k_window = zeros(Float64, x, y)
    for i = 1:x:N_x*x-(x-1), j = 1:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        k[div(i, x)+1, div(j, y)+1] = geomean(k_window)
    end
    return k
end

function harmonic_average(k_old::Array, x::Int64, y::Int64)
    N_x = size(k_old, 1) ÷ x
    N_y = size(k_old, 2) ÷ y
    k = zeros(N_x, N_y)
    k_window = zeros(Float64, x, y)
    for i = 1:x:N_x*x-(x-1), j = 1:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        k[div(i, x)+1, div(j, y)+1] = harmmean(k_window)
    end
    return k
end

function power_average(k_old::Array, x::Int64, y::Int64)
    N_x = size(k_old, 1) ÷ x
    N_y = size(k_old, 2) ÷ y
    k = zeros(N_x, N_y)
    k_window = zeros(Float64, x, y)
    for i = 1:x:N_x*x-(x-1), j = 1:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        k[div(i, x)+1, div(j, y)+1] = exp(sum(log.(k_window)) / (x * y))
    end
    return k
end

function simple_average(k_old::Array, x::Int64, y::Int64)
    N_x = size(k_old, 1) ÷ x
    N_y = size(k_old, 2) ÷ y
    k = zeros(N_x, N_y)
    k_window = zeros(Float64, x, y)
    for i = 1:x:N_x*x-(x-1), j = 1:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        k[div(i, x)+1, div(j, y)+1] = mean(k_window)
    end
    return k
end
