using Statistics, Plots
include("gradient.jl")
include("p_upsc.jl")
include("fvm_1d.jl")

N = 10
L_x = 10 #расстояние (м)
P_begin = 5*10^6; P_end = 10^6
h_x = L_x / N
dx = 10 / N
k1_old = rand(Float64, N + 2, 1)
k2_old = upscaling(k1_old, N + 2)
mu = 10^(-3)
h = 0.01
k1 = zeros(N, 1)
k2 = zeros(2 * N - 1, 1)
u_val_1 = zeros(N, 1)
u_val_2 = zeros(2 * N - 1, 1)

p1 = FVM(N, P_begin, P_end, h_x, k1_old, mu, h)
p2 = FVM(2*N-1, P_begin, P_end, h_x, k2_old, mu, h)

for i = 1:N
    k1[i] = k1_old[i+1]
end

for i = 1:2*N-1
    k2[i] = k2_old[i+1]
end

for i = 1:N
    u_val_1[i] = - (k1[i] / mu) * grad1D(p1, dx)[i] * 10^(-12)
end

for i = 1:2*N-1
    u_val_2[i] = - (k2[i] / mu) * grad1D(p2, dx)[i] * 10^(-12)
end

sr_arifm_1 = sum(u_val_1) / N
sr_arifm_2 = sum(u_val_2) / (2*N-1)
otklon_1 = std(u_val_1)
otklon_2 = std(u_val_2)
sr = mean(u_val_1)
println("среднее арифметическое первое = ", sr_arifm_1)
println("среднее арифметическое второе = ", sr_arifm_2)
println("ср.квадратичн. отклон. первое = ", otklon_1)
println("ср.квадратичн. отклон. второе = ", otklon_2)

x = 1:N+1
x1 = 1:N
x2 = 1:2*N-1
x3 = 1:2*N

plotly()
#Plots.scalefontsizes(1.2)

fig1 = plot(x1, k1, width=3)
xlabel!("X, cm")
ylabel!("k, D")

fig2 = plot(x, p1, linewidth=3)
xlabel!("X, cm")
ylabel!("P, Pa")

fig3 = plot(x1, u_val_1, linewidth=3)
xlabel!("X, cm")
ylabel!("U, m/s")

fig4 = plot(x2/2, k2, linewidth=3)
xlabel!("X, cm")
ylabel!("k, D")

fig5 = plot(x3/2, p2, linewidth=3)
xlabel!("X, cm")
ylabel!("P, Pa")

fig6 = plot(x2/2, u_val_2, linewidth=3)
xlabel!("X, cm")
ylabel!("U, m/s")

display(plot(fig1, fig2, fig3, layout=(3,1)))
display(plot(fig4, fig5, fig6, layout=(3,1)))
