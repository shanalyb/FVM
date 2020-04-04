using PyCall, Statistics
import PyPlot: quiver
using PyPlot
pygui(true)
include("gradient.jl")
include("concentration_2d.jl")
include("k_generator.jl")
include("upscaling.jl")
include("pressure.jl")
include("velocity.jl")
N_x = 50
N_y = 50
L_x = 10 #Расстояние (м)
L_y = 10 #Расстояние (м)
P_left = 4.0
P_right = 3.9
P_up = 4.0
P_down = 4.0
h_x = L_x / N_x
h_y = L_y / N_y
dx = 10 / N_x
dy = 10 / N_y
#k = rand(Float64, N_x + 2, N_y + 2)
mu = 10^(-3)
h = 0.01
#k = 0.5
x = 6
y = 6
#Генерация проницаемости
#--------------------------------
range_1 = 3
alpha_1 = 1.0
beta_1 = 5.0
tetha_1 = -pi / 4
S = 0.5
k =
    permeability_generator(
        N_x + 3,
        N_y + 3,
        range_1,
        alpha_1,
        beta_1,
        tetha_1,
        S,
    ) * 10
k_harmonic = harmonic_average(k, x, y)
k_geometric = geometric_average(k, x, y)
k_power = power_average(k, x, y)
#--------------------------------

p = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, h, k)
p_harmonic = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, h, k_harmonic)
p_geometric = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, h, k_geometric)
p_power = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, h, k_power)

#Скорости в СИ
u_x = velocity_x(k, mu, p, h_x) * 10^(-6)
u_y = velocity_y(k, mu, p, h_y) * 10^(-6)
u_summ = sqrt.(u_x .^ 2 + u_y .^ 2)
koef = u_y ./ u_x

u_x_harmonic = velocity_x(k_harmonic, mu, p_harmonic, h_x) * 10^(-6)
u_y_harmonic = velocity_y(k_harmonic, mu, p_harmonic, h_y) * 10^(-6)
u_harmonic = sqrt.(u_x_harmonic .^ 2 + u_y_harmonic .^ 2)
koef_harmonic = u_y_harmonic ./ u_x_harmonic

u_x_geometric = velocity_x(k_geometric, mu, p_geometric, h_x) * 10^(-6)
u_y_geometric = velocity_y(k_geometric, mu, p_geometric, h_y) * 10^(-6)
u_geometric = sqrt.(u_x_geometric .^ 2 + u_y_geometric .^ 2)
koef_geometric = u_y_geometric ./ u_x_geometric

u_x_power = velocity_x(k_power, mu, p_power, h_x) * 10^(-6)
u_y_power = velocity_y(k_power, mu, p_power, h_y) * 10^(-6)
u_power = sqrt.(u_x_power .^ 2 + u_y_power .^ 2)
koef_power = u_y_power ./ u_x_power
#Концентрация примеси
#--------------------------------
c_0 = zeros(size(u_x, 1) + 1, size(u_y, 2) + 1)
for i = 1:size(u_x, 1) + 1
    c_0[i, 1] = 1
end
c_0_harmonic = zeros(size(u_x_harmonic, 1) + 1, size(u_y_harmonic, 2) + 1)
for i = 1:size(u_x_harmonic, 1) + 1
    c_0_harmonic[i, 1] = 1
end
c_0_geometric = zeros(size(u_x_geometric, 1) + 1, size(u_y_geometric, 2) + 1)
for i = 1:size(u_x_geometric, 1) + 1
    c_0_geometric[i, 1] = 1
end
c_0_power = zeros(size(u_x_power, 1) + 1, size(u_y_power, 2) + 1)
for i = 1:size(u_x_power, 1) + 1
    c_0_power[i, 1] = 1
end

t = 170
c = concentration_2d(c_0, h_x, h_y, u_x, u_y, t)
c_harmonic = concentration_2d(c_0_harmonic, h_x, h_y, u_x_harmonic, u_y_harmonic, t)
c_geometric = concentration_2d(c_0_geometric, h_x, h_y, u_x_geometric, u_y_geometric, t)
c_power = concentration_2d(c_0_power, h_x, h_y, u_x_power, u_y_power, t)
#--------------------------------
