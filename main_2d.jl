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
include("tensor_visualization.jl")
N_x = 50
N_y = 50
L_x = 10 #Расстояние (м)
L_y = 10 #Расстояние (м)
P_left = 4.0
P_right = 2.0
P_up = 0.0
P_down = 0.0
h_x = L_x / N_x
h_y = L_y / N_y
dx = 10 / N_x
dy = 10 / N_y
#k = rand(Float64, N_x + 2, N_y + 2)
mu = 10^(-3)
h = 0.01
#k_1 = 0.5
x = 3
y = 3
#Генерация проницаемости
#--------------------------------
range_1 = 3
alpha_1 = 1.0
beta_1 = 1.0
tetha_1 = 0.0
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
    )/10
k_harmonic = harmonic_average(k, x, y)
k_geometric = geometric_average(k, x, y)
k_power = power_average(k, x, y)
k_generalized = generalized_average(k, x, y)
#--------------------------------

p = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, k)
p_harmonic = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_harmonic)
p_geometric = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_geometric)
p_power = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_power)
p_generalized = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_generalized)


p = p[setdiff(1:end, (1, N_x + 2)), setdiff(1:end, (1, N_y + 2))]
#Скорости в СИ
u_x = velocity_x(k, mu, p, h_x) * 10^(-6)
u_y = velocity_y(k, mu, p, h_y) * 10^(-6)
u_summ = sqrt.(u_x .^ 2 + u_y .^ 2)
koef = u_y ./ u_x

u_x_harmonic = velocity_x(k_harmonic, mu, p_harmonic, h_x*x) * 10^(-6)
u_y_harmonic = velocity_y(k_harmonic, mu, p_harmonic, h_y*y) * 10^(-6)
u_harmonic = sqrt.(u_x_harmonic .^ 2 + u_y_harmonic .^ 2)
koef_harmonic = u_y_harmonic ./ u_x_harmonic

u_x_geometric = velocity_x(k_geometric, mu, p_geometric, h_x*x) * 10^(-6)
u_y_geometric = velocity_y(k_geometric, mu, p_geometric, h_y*y) * 10^(-6)
u_geometric = sqrt.(u_x_geometric .^ 2 + u_y_geometric .^ 2)
koef_geometric = u_y_geometric ./ u_x_geometric

u_x_power = velocity_x(k_power, mu, p_power, h_x*x) * 10^(-6)
u_y_power = velocity_y(k_power, mu, p_power, h_y*y) * 10^(-6)
u_power = sqrt.(u_x_power .^ 2 + u_y_power .^ 2)
koef_power = u_y_power ./ u_x_power

u_x_generalized = velocity_x(k_generalized, mu, p_generalized, h_x*x) * 10^(-6)
u_y_generalized = velocity_y(k_generalized, mu, p_generalized, h_y*y) * 10^(-6)
u_generalized = sqrt.(u_x_generalized .^ 2 + u_y_generalized .^ 2)
koef_generalized = u_y_generalized ./ u_x_generalized
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
c_0_generalized = zeros(size(u_x_generalized, 1) + 1, size(u_y_generalized, 2) + 1)
for i = 1:size(u_x_generalized, 1) + 1
    c_0_generalized[i, 1] = 1
end

t_main = 200
c = concentration_2d(c_0, h_x, h_y, u_x, u_y, t_main)
c_harmonic = concentration_2d(c_0_harmonic, h_x*x, h_y*y, u_x_harmonic, u_y_harmonic, t_main)
c_geometric = concentration_2d(c_0_geometric, h_x*x, h_y*y, u_x_geometric, u_y_geometric, t_main)
c_power = concentration_2d(c_0_power, h_x*x, h_y*y, u_x_power, u_y_power, t_main)
c_generalized = concentration_2d(c_0_generalized, h_x*x, h_y*y, u_x_generalized, u_y_generalized, t_main)
#--------------------------------


#u_x_transmissibility = transmissibility_average(k::Array, x::Int64, y::Int64, mu::Float64)[1,1] * 10^(-6)
#u_y_transmissibility = transmissibility_average(k::Array, x::Int64, y::Int64, mu::Float64)[2,1] * 10^(-6)

k_x = transmissibility_average(k, x, y, 1.0)[1,1]
k_y = transmissibility_average(k, x, y, 1.0)[2,2]
k_transmissibility = 0.5*(k_x + k_y)

p_x_transmissibility = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_x)
p_y_transmissibility = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_y)
p_transmissibility = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_transmissibility)

u_x_transmissibility = velocity_x(k_x, mu, p_x_transmissibility, h_x*x) * 10^(-6)
u_y_transmissibility = velocity_y(k_y, mu, p_y_transmissibility, h_y*y) * 10^(-6)
u_transmissibility = sqrt.(u_x_transmissibility .^ 2 + u_y_transmissibility .^ 2)
koef_transmissibility = u_y_transmissibility ./ u_x_transmissibility

c_0_transmissibility = zeros(size(u_x_transmissibility, 1) + 1, size(u_y_transmissibility, 2) + 1)
for i = 1:size(u_x_transmissibility, 1) + 1
    c_0_transmissibility[i, 1] = 1
end

c_transmissibility = concentration_2d(c_0_transmissibility, h_x*x, h_y*y, u_x_transmissibility, u_y_transmissibility, t_main)


k_new = zeros(3, 3)
k_new[1,1]=k_transmissibility[1,1][1]
k_new[2,2]=k_transmissibility[2,2][1]
#fig = plt.figure()
#fig.suptitle("Original and power average", fontsize=14)

#ax = plt.subplot(1,1,1)
#ax.set_title("pressure and velocity")
#im_2 = ax.imshow(p_generalized, extent = [0, size(p_generalized, 1)-1, 0, size(p_generalized, 2)-1], interpolation = "gaussian")
#cbar = colorbar(im_2)
#cbar.set_label("P, MPa", fontsize = 14)
#quiver(u_x_transmissibility, u_y_transmissibility)


fig = plt.figure()
fig.suptitle("Концентрация примеси", fontsize=14)
#Концентрация
ax = plt.subplot(1,1,1)
im_3 =ax.imshow(c, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.ax.tick_params(labelsize=20)

plt.tick_params(labelsize=20)
plt.subplots_adjust(wspace=0, hspace=0.3)
plt.show()

#draw_tensor(k_new, "example.png")
