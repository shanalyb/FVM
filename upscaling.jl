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

function generalized_average(k_old::Array, x::Int64, y::Int64)
    N_x = size(k_old, 1) ÷ x
    N_y = size(k_old, 2) ÷ y
    k = zeros(N_x, N_y)
    k_window = zeros(Float64, x, y)
    for i = 1:x:N_x*x-(x-1), j = 1:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        k[div(i, x)+1, div(j, y)+1] = genmean(k_window, 2)
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

function transmissibility_average(k_old::Array, x::Int64, y::Int64, p::Array, mu::Float64)
    N_x = size(p, 1) ÷ x
    N_y = size(p, 2) ÷ y
    k_x = zeros(N_x-1, N_y-1)
    k_y = zeros(N_x-1, N_y-1)
    grad_p_x = grad_p_y = zeros(N_x-1, N_y-1)
    k_window = zeros(Float64, x, y)
    p_window = zeros(Float64, x+1, y+1)

    for i = 2:x:N_x*x-(x-1), j = 2:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        for m = 1:x+1, k = 1:y+1
            p_window[m, k] = p[i+(m-1)-1, j+(k-1)-1]
        end
        grad_p_x[div(i-1, x)+1, div(j-1, y)+1] = mean(grad2D_x(p_window, h_x))
        k_x[div(i-1, x)+1, div(j-1, y)+1] = mean(k_window .* grad2D_x(p_window, h_x)) / mean(grad2D_x(p_window, h_x))
    end

    for i = 2:x:N_x*x-(x-1), j = 2:y:N_y*y-(y-1)
        for m = 1:x, k = 1:y
            k_window[m, k] = k_old[i+(m-1), j+(k-1)]
        end
        for m = 1:x+1, k = 1:y+1
            p_window[m, k] = p[i+(m-1)-1, j+(k-1)-1]
        end
        grad_p_y[div(i-1, x)+1, div(j-1, y)+1] = mean(grad2D_y(p_window, h_y))
        k_y[div(i-1, x)+1, div(j-1, y)+1] = mean(k_window .* grad2D_y(p_window, h_y)) / mean(grad2D_y(p_window, h_y))
    end

    u_x = - (1 / mu) * (k_x .* grad_p_x)
    u_y = - (1 / mu) * (k_y .* grad_p_y)

    u = Array{Array,2}(undef, 2, 1)
    u[1, 1] = u_x
    u[2, 1] = u_y
    return u
end
