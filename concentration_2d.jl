function concentration_2d(c_0, N_x, N_y, h_x, h_y, u_x, u_y, delta_t)
    dif_c_x = grad2D_x(c_0, h_x)
    dif_c_y = grad2D_y(c_0, h_y)
    c = zeros(N_x+1, N_y+1)

    dif_c_x1 = zeros(N_x+1, N_y+1)
    for i = 1:N_x+1, j = 1:N_y
        dif_c_x1[i,j] = dif_c_x[i,j]
    end

    dif_c_y1 = zeros(N_x+1, N_y+1)
    for i = 1:N_x, j = 1:N_y+1
        dif_c_y1[i,j] = dif_c_y[i,j]
    end

    u_x1 = zeros(N_x+1, N_y+1)
    for i = 1:N_x+1, j = 1:N_y
        u_x1[i,j] = u_x[i,j]
    end

    u_y1 = zeros(N_x+1, N_y+1)
    for i = 1:N_x, j = 1:N_y+1
        u_y1[i,j] = u_y[i,j]
    end

    c = c_0 - (u_y1 * dif_c_x1 + u_x1 * dif_c_y1) * delta_t
end
