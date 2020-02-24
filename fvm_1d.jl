function FVM(N, P_begin, P_end, h_x, k, mu, h)
a_old = zeros(N + 1, 1)
a = zeros(N + 1, N + 1)
b = zeros(N + 1, 1)

for i = 2:N+2
    a_old[i-1] = (h / (k[i] - k[i-1])) * (log(k[i] / k[i-1]))
end

for k = 1:N-1
    a[k, k] = -a_old[k]
    a[k, k+1] = a_old[k+1] + a_old[k]
    a[k, k+2] = -a_old[k+1]
end
a[N, 1] = 1
a[N+1, N+1] = 1

a_new = inv(a)

b[N] = P_begin
b[N+1] = P_end

p = a_new * b
end
