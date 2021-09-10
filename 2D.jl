using PyCall, Statistics
import PyPlot: quiver
using PyPlot
pygui(true)
include("gradient.jl")
include("k_generator.jl")
include("pressure.jl")
include("velocity.jl")
N_x = 25
N_y = 25
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
#k = 0.5
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

p = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x, h_y, mu, k)

u_x = velocity_x(k, mu, p, h_x) * 10^(-6)
u_y = velocity_y(k, mu, p, h_y) * 10^(-6)

fig = plt.figure()
#fig.suptitle("Original and power average", fontsize=14)

x = range(0, stop=250, length=26)
y = range(250, stop=0, length=26)


fig = plt.figure()
ax = plt.subplot(1,1,1)
#ax.set_title("pressure and velocity")
im_2 = ax.imshow(p, extent = [0, 250, 0, 250], interpolation = "gaussian")
cbar = colorbar(im_2)
cbar.set_label(raw"$P, MPa$", fontsize = 16)
cbar.ax.tick_params(labelsize=16)
quiver(x, y, u_x, u_y)

plt.tick_params(labelsize=16)
plt.subplots_adjust(wspace=0, hspace=0.3)
plt.show()


fig = plt.figure()
ax = plt.subplot(1,1,1)
#ax.set_title("pressure and velocity")
ax.set_title("Поле проницаемости", y=1.1)
im_1 = ax.imshow(k, extent = [0, 250, 0, 250])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 24)
cbar.ax.tick_params(labelsize=24)

plt.tick_params(labelsize=24)
plt.subplots_adjust(wspace=0, hspace=0.3)
plt.show()
