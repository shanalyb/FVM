function pressure_homogeneous_medium(
    P_left::Float64,
    P_right::Float64,
    P_up::Float64,
    P_down::Float64,
    h_x::Float64,
    h_y::Float64,
    mu::Float64,
    h::Float64,
    k::Array
)
    N_x = size(k, 1)-1
    N_y = size(k, 2)-1
    a_old = zeros(N_x + 1, N_y + 1)
    b_old = zeros(N_x + 1, N_y + 1)
    A = zeros((N_x + 1) * (N_y + 1), (N_x + 1) * (N_y + 1))
    b = zeros((N_x + 1) * (N_y + 1), 1)
    p = zeros(N_x + 1, N_y + 1)
    p_min = zeros(N_y + 1, 1)
    p_max = zeros(N_y + 1, 1)
    p_ex = zeros(1, (N_x + 1) * (N_y + 1))

    for i = 1:N_x+1, j = 1:N_x+1
        a_old[i, j] = (k / h_x)
    end

    for i = 1:N_y+1, j = 1:N_y+1
        b_old[i, j] = (k / h_y)
    end

    for m = 1:(N_y-1)
        A[m, m+1] = 1
        A[N_y-1+m, N_x * (N_y + 1)+m+1] = 1
        b[m, 1] = P_left
        b[N_y-1+m, 1] = P_right
    end

    for m = 0:N_x
        A[2*(N_y-1)+(m+1), m*(N_y+1)+1] = 1
        A[2 * (N_y - 1)+(N_x + 1)+(m + 1), m*(N_y+1)+(N_y+1)] = 1
        #b[2*(N_y-1)+(m+1), 1] = P_up
        #b[2*(N_y-1)+(N_x+1)+(m+1), 1] = P_down
        A[2*(N_y-1)+(m+1), m * (N_y + 1)+1+1] = -1
        A[2 * (N_y - 1)+(N_x + 1)+(m + 1), m*(N_y+1)+(N_y+1)-1] = -1
    end

    m = 2 * (N_x + N_y)
    for i = 1:N_x-1, j = 1:N_y-1
        m = m + 1
        A[m, (i-1)*(N_y+1)+(j+1)] = -a_old[i, j]
        A[m, i*(N_y+1)+(j+1)] =
            a_old[i, j] + a_old[i+1, j] + b_old[i, j] + b_old[i, j+1]
        A[m, (i+1)*(N_y+1)+(j+1)] = -a_old[i+1, j]
        A[m, i*(N_y+1)+j] = -b_old[i, j]
        A[m, i*(N_y+1)+(j+2)] = -b_old[i, j+1]
    end

    A_new = inv(A)

    p_old = A_new * b

    n = 0
    for j = 1:N_y+1, i = 1:N_x+1
        n = n + 1
        p[i, j] = p_old[n]
    end
    return p
end

function pressure_inhomogeneous_medium(
    P_left::Float64,
    P_right::Float64,
    P_up::Float64,
    P_down::Float64,
    h_x::Float64,
    h_y::Float64,
    mu::Float64,
    h::Float64,
    k::Array
)
    N_x = size(k, 1)-1
    N_y = size(k, 2)-1
    a_old = zeros(N_x + 1, N_y + 1)
    b_old = zeros(N_x + 1, N_y + 1)
    A = zeros((N_x + 1) * (N_y + 1), (N_x + 1) * (N_y + 1))
    b = zeros((N_x + 1) * (N_y + 1), 1)
    p = zeros(N_x + 1, N_y + 1)
    p_min = zeros(N_y + 1, 1)
    p_max = zeros(N_y + 1, 1)
    p_ex = zeros(1, (N_x + 1) * (N_y + 1))

    #Неоднородная среда
    #--------------------------------
    for i = 2:size(k, 1), j = 1:size(k, 2)-1
        b_old[i-1, j] =
            ((h_y / (k[i, j] - k[i-1, j])) * (log(k[i, j] / k[i-1, j])))^(-1)
    end

    for i = 1:size(k, 1)-1, j = 2:size(k, 2)
        a_old[i, j-1] =
            ((h_x / (k[i, j] - k[i, j-1])) * (log(k[i, j] / k[i, j-1])))^(-1)
    end
    #--------------------------------

    for m = 1:(N_y-1)
        A[m, m+1] = 1
        A[N_y-1+m, N_x * (N_y + 1)+m+1] = 1
        b[m, 1] = P_left
        b[N_y-1+m, 1] = P_right
    end

    for m = 0:N_x
        A[2*(N_y-1)+(m+1), m*(N_y+1)+1] = 1
        A[2 * (N_y - 1)+(N_x + 1)+(m + 1), m*(N_y+1)+(N_y+1)] = 1
        #b[2*(N_y-1)+(m+1), 1] = P_up
        #b[2*(N_y-1)+(N_x+1)+(m+1), 1] = P_down
        A[2*(N_y-1)+(m+1), m * (N_y + 1)+1+1] = -1
        A[2 * (N_y - 1)+(N_x + 1)+(m + 1), m*(N_y+1)+(N_y+1)-1] = -1
    end

    m = 2 * (N_x + N_y)
    for i = 1:N_x-1, j = 1:N_y-1
        m = m + 1
        A[m, (i-1)*(N_y+1)+(j+1)] = -a_old[i, j]
        A[m, i*(N_y+1)+(j+1)] =
            a_old[i, j] + a_old[i+1, j] + b_old[i, j] + b_old[i, j+1]
        A[m, (i+1)*(N_y+1)+(j+1)] = -a_old[i+1, j]
        A[m, i*(N_y+1)+j] = -b_old[i, j]
        A[m, i*(N_y+1)+(j+2)] = -b_old[i, j+1]
    end

    A_new = inv(A)

    p_old = A_new * b
    p = reshape(p_old, (N_x+1, N_y+1))
    #n = 0
    #for j = 1:N_y+1, i = 1:N_x+1
    #    n = n + 1
    #    p[i, j] = p_old[n]
    #end
    return p
end
