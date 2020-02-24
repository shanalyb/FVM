#function grad1D(P::Array, k::Float64)
#    Nj = P.domain.dims[1] + 1
#    Ni = P.domain.dims[2] + 1
#    for j = 2:Nj, i = 2:Ni
#        U.value[i-1, j] = k * (-2 * P.value[i, j] + P.value[i-1, j] +
#                           P.value[i+1, j])
#    end
#    U
#end

#function grad1D_first(P::Array, k::Float64)
#    Ni = P.domain.dims[1] + 1
#    for i = 2:Ni
#        U.value[i-1] = k * (-2 * P.value[i] + P.value[i-1] +
#                           P.value[i+1])
#    end
#    U
#end

function grad1D(P::Array, h_x::Float64)
    Ni = size(P, 1)
    U = zeros(Ni - 1, 1)
    for i = 2:Ni
        U[i-1] = (P[i] - P[i-1]) / h_x
    end
    U
end

function grad2D_x(P::Array, h_x::Float64)
    N_x = size(P, 1); N_y = size(P, 2) #1-по вертикали, 2-по горизонтали
    U_x = zeros(N_x, N_y-1)
    for i = 1:N_x, j = 1:N_y-1
        U_x[i, j] = (P[i, j+1] - P[i, j]) / h_x
    end
    U_x
end

function grad2D_y(P::Array, h_y::Float64)
    N_x = size(P, 1); N_y = size(P, 2)
    U_y = zeros(N_x-1, N_y)
    for j = 1:N_y, i = 1:N_x-1
        U_y[i, j] = (P[i+1, j] - P[i, j]) / h_y
    end
    U_y
end
