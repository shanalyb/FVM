function FVM(N_x, N_y, P_left, P_right, P_up, P_down, h_x, k_a, k_b, mu, h_x, h_y)
a_old = zeros(N + 1, 1)
A = zeros(6 * (N_x + 1), N_y + 1)
b = zeros((N_x + 1)*(N_y + 1), 1)

for i = 2:N_x+2, j = 1:N_x+1
    a_old[i-1, j] = (h_y / (k[i, j] - k[i-1, j])) * (log(k[i, j] / k[i-1, j]))
end

for i = 1:N_y+1, j = 2:N_y+2
    b_old[i, j-1] = (h_x / (k[i, j] - k[i, j-1])) * (log(k[i, j] / k[i, j-1]))
end

for i = 1:(N_x-1)*(N_y+1)-N_x
    A[i, i+1] = -a_old[i, j]
    A[i, i+N_x+1] = b_old[i, j]
    A[i, i+N_x+2] = a_old[i, j] + a_old[i+1, j] + b_old[i, j] + b_old[i, j+1]
    A[i, i+N_x+3] = -b_old[i, j+1]
    A[i, i+2*N_x+2] = a_old[i+1, j]
end

for j = 1:N_y
    A[(N_x-1)*(N_y+1)-N_x+1, j] = 1
end

for i = (N_x-1)*(N_y+1)-N_x+1:(N_x-1)*(N_y+1), j = N_y+1
    A[i, j] = 1

a_new = inv(a)

b[N] = X_begin
b[N+1] = X_end

p = a_new * b
end
