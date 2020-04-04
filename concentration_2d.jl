function concentration_2d(c_0::Array, h_x::Float64, h_y::Float64, u_x::Array, u_y::Array, t::Int64)
    N_x = size(c_0, 1) - 1
    N_y = size(c_0, 2) - 1
    for i = 1:t
        c_x = zeros(N_x + 1, N_y)
        c_y = zeros(N_x, N_y + 1)
        dif_c_x = grad2D_x(c_0, h_x)
        dif_c_y = grad2D_y(c_0, h_y)

        dif_c_x1 = zeros(N_x + 1, N_y + 1)
        for i = 1:N_x, j = 1:N_y
            dif_c_x1[i, j+1] = dif_c_x[i, j]
        end

        dif_c_y1 = zeros(N_x + 1, N_y + 1)
        for i = 1:N_x, j = 1:N_y
            dif_c_y1[i, j] = dif_c_y[i, j]
        end

        for i = 1:N_x, j = 1:N_y
            c_x[i, j] = u_x[i, j] * dif_c_x1[i, j]
        end
        for i = 1:N_x, j = 1:N_y
            c_y[i, j] = u_y[i, j] * dif_c_y1[i, j]
        end

        c_x1 = zeros(N_x + 1, N_y + 1)
        for i = 1:N_x+1, j = 1:N_y
            c_x1[i, j] = c_x[i, j]
        end

        c_y1 = zeros(N_x + 1, N_y + 1)
        for i = 1:N_x, j = 1:N_y+1
            c_y1[i, j] = c_y[i, j]
        end

        c_0 = c_0 - (c_x1 + c_y1) * t
    end

    c = zeros(N_x, N_y)
    for i = 1:N_x, j = 1:N_y
        c[i, j] = c_0[i, j]
    end
    return c
end
