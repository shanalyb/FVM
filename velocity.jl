function velocity_x(k, mu, p, h_x)
    N_x = size(p, 1) - 1
    N_y = size(p, 2) - 1
    k_0 = zeros(N_x, N_y)
    for i = 1:N_x, j = 1:N_y
        k_0[i, j] = k[i, j]
    end
    u_x = -(k_0 / mu) .* grad2D_x(p, h_x)
end

function velocity_y(k, mu, p, h_y)
    N_x = size(p, 1) - 1
    N_y = size(p, 2) - 1
    k_0 = zeros(N_x, N_y)
    for i = 1:N_x, j = 1:N_y
        k_0[i, j] = k[i, j]
    end
    u_y = -(k_0 / mu) .* grad2D_y(p, h_y)
end
