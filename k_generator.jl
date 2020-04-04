function permeability_generator(
    N_x::Int64,
    N_y::Int64,
    range::Int64,
    alpha::Float64,
    beta::Float64,
    tetha::Float64,
    S::Float64,
)
    a = rand(Float64, N_x + 1, N_y + 1)
    a_new = reshape(a, ((N_x + 1) * (N_y + 1), 1))
    A = zeros((N_x + 1) * (N_y + 1), (N_x + 1) * (N_y + 1))
    for i = 0:N_x, j = 0:N_y, l = 0:N_x, m = 0:N_y
        k = i * (N_y + 1) + (j + 1)
        n = l * (N_y + 1) + (m + 1)
        alpha_new = alpha*cos(tetha) + beta*sin(tetha)
        beta_new = -alpha*sin(tetha) + beta*cos(tetha)
        x = alpha * (cos(tetha) * (i - l) + sin(tetha) * (j - m))
        y = beta * (-sin(tetha) * (i - l) + cos(tetha) * (j - m))
        r = sqrt(x^2 + y^2)
        if r <= range
            A[k, n] = ((range - r) / range)^S
        else
            A[k, n] = 0
        end
    end
    nu = A * a_new
    nu_matrix = reshape(nu, (N_x + 1, N_y + 1))
    k = nu_matrix[setdiff(1:end, (1, N_x + 1)), setdiff(1:end, (1, N_y + 1))]
    return k
end
