using PyCall, DataFrames
import PyPlot: quiver
using PyPlot
pygui(true)
include("main_2d.jl")
include("upscaling.jl")
k_simple = simple_average(k, x, y)
p_simple = simple_average(p, x, y)
u_simple = simple_average(u_summ, x, y)[setdiff(1:end, (1,end-1,end)), setdiff(1:end, (1,end-1,end))]
u_x_simple = simple_average(u_x, x, y)[setdiff(1:end, (1,end-1,end)), setdiff(1:end, (1,end-1,end))]
u_y_simple = simple_average(u_y, x, y)[setdiff(1:end, (1,end-1,end)), setdiff(1:end, (1,end-1,end))]


p_harmonic_1 = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_harmonic[setdiff(1:end, 1), setdiff(1:end, 1)])
p_geometric_1 = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_geometric[setdiff(1:end, 1), setdiff(1:end, 1)])
p_power_1 = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_power[setdiff(1:end, 1), setdiff(1:end, 1)])
p_generalized_1 = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_generalized[setdiff(1:end, 1), setdiff(1:end, 1)])
p_transmissibility_1 = pressure_inhomogeneous_medium(P_left, P_right, P_up, P_down, h_x*x, h_y*y, mu, k_transmissibility)



k_1 = k[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]
#p_1 = p[setdiff(1:end, 1), setdiff(1:end, 1)]
#u_summ_1 = u_summ[setdiff(1:end, 1), setdiff(1:end, 1)]

k_simple_1 = simple_average(k_1, x, y)
#p_simple_1 = simple_average(p_1, x, y)
#u_simple_1 = simple_average(u_summ_1, x, y)

p_simple_1 = p_simple[setdiff(1:end, (1,end-2,end-1,end)), setdiff(1:end, (1,end-2,end-1,end))]
#p_simple = p_simple[setdiff(1:end, end), setdiff(1:end, end)]
#k_simple_1 = k_simple[setdiff(1:end ,1), setdiff(1:end ,1)]
u_simple_1 = u_simple[setdiff(1:end, end), setdiff(1:end, end)]

function equalize_the_sizes(p, p_small)
    p_new = zeros(size(p_small, 1), size(p_small, 2))
    for i = 1:size(p_small, 1), j = 1:size(p_small, 2)
        p_new[i, j] = p[i+1, j+1]
    end
    return p_new
end


p_harmonic_1 = p_harmonic_1[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]
p_geometric_1 = p_geometric_1[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]
p_power_1 = p_power_1[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]
p_generalized_1 = p_generalized_1[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]
p_transmissibility_1 = p_transmissibility_1[setdiff(1:end, (1,end)), setdiff(1:end, (1,end))]

u_x_harmonic = velocity_x(k_harmonic, mu, p_harmonic_1, h_x*x) * 10^(-6)
u_y_harmonic = velocity_y(k_harmonic, mu, p_harmonic_1, h_y*y) * 10^(-6)
u_harmonic = sqrt.(u_x_harmonic .^ 2 + u_y_harmonic .^ 2)

u_x_geometric = velocity_x(k_geometric, mu, p_geometric_1, h_x*x) * 10^(-6)
u_y_geometric = velocity_y(k_geometric, mu, p_geometric_1, h_y*y) * 10^(-6)
u_geometric = sqrt.(u_x_geometric .^ 2 + u_y_geometric .^ 2)

u_x_power = velocity_x(k_power, mu, p_power_1, h_x*x) * 10^(-6)
u_y_power = velocity_y(k_power, mu, p_power_1, h_y*y) * 10^(-6)
u_power = sqrt.(u_x_power .^ 2 + u_y_power .^ 2)

u_x_generalized = velocity_x(k_generalized, mu, p_generalized_1, h_x*x) * 10^(-6)
u_y_generalized = velocity_y(k_generalized, mu, p_generalized_1, h_y*y) * 10^(-6)
u_generalized = sqrt.(u_x_generalized .^ 2 + u_y_generalized .^ 2)

u_x_transmissibility = velocity_x(k_x, mu, p_transmissibility_1, h_x*x) * 10^(-6)
u_y_transmissibility = velocity_y(k_y, mu, p_transmissibility_1, h_y*y) * 10^(-6)
u_transmissibility = sqrt.(u_x_transmissibility .^ 2 + u_y_transmissibility .^ 2)


p_harmonic_1 = p_harmonic_1[setdiff(1:end, (end-1, end)), setdiff(1:end, (end-1, end))]
p_geometric_1 = p_geometric_1[setdiff(1:end, (end-1, end)), setdiff(1:end, (end-1, end))]
p_power_1 = p_power_1[setdiff(1:end, (end-1, end)), setdiff(1:end, (end-1, end))]
p_generalized_1 = p_generalized_1[setdiff(1:end, (end-1, end)), setdiff(1:end, (end-1, end))]
p_transmissibility_1 = p_transmissibility_1[setdiff(1:end, (end-1, end)), setdiff(1:end, (end-1, end))]
#u_geometric = equalize_the_sizes(u_geometric, u_simple)
#u_harmonic = equalize_the_sizes(u_harmonic, u_simple)
#u_power = equalize_the_sizes(u_power, u_simple)
#u_generalized = equalize_the_sizes(u_generalized, u_simple)

k_difference_harmonic = abs.(k_simple - k_harmonic)
k_difference_geometric = abs.(k_simple - k_geometric)
k_difference_power = abs.(k_simple - k_power)
k_difference_generalized = abs.(k_simple - k_generalized)
k_difference_transmissibility = abs.(k_simple_1 - k_transmissibility)

p_difference_harmonic = abs.(p_simple_1 - p_harmonic_1)
p_difference_geometric = abs.(p_simple_1 - p_geometric_1)
p_difference_power = abs.(p_simple_1 - p_power_1)
p_difference_generalized = abs.(p_simple_1 - p_generalized_1)
p_difference_transmissibility = abs.(p_simple_1 - p_transmissibility_1)

u_difference_harmonic = abs.(u_simple - u_harmonic)
u_difference_geometric = abs.(u_simple - u_geometric)
u_difference_power = abs.(u_simple - u_power)
u_difference_generalized = abs.(u_simple - u_generalized)
u_difference_transmissibility = abs.(u_simple - u_transmissibility)

#Сравнение отклонений
k_difference_harmonic_max = maximum(k_difference_harmonic)
k_difference_geometric_max = maximum(k_difference_geometric)
k_difference_power_max = maximum(k_difference_power)
k_difference_generalized_max = maximum(k_difference_generalized)
k_difference_transmissibility_max = maximum(k_difference_transmissibility)

p_difference_harmonic_max = maximum(p_difference_harmonic)
p_difference_geometric_max = maximum(p_difference_geometric)
p_difference_power_max = maximum(p_difference_power)
p_difference_generalized_max = maximum(p_difference_generalized)
p_difference_transmissibility_max = maximum(p_difference_transmissibility)

u_difference_harmonic_max = maximum(u_difference_harmonic)
u_difference_geometric_max = maximum(u_difference_geometric)
u_difference_power_max = maximum(u_difference_power)
u_difference_generalized_max = maximum(u_difference_generalized)
u_difference_transmissibility_max = maximum(u_difference_transmissibility)

df = DataFrame(parameters = ["k", "p", "u"], harmonic = [k_difference_harmonic_max, p_difference_harmonic_max, u_difference_harmonic_max],
                                            geometric = [k_difference_geometric_max, p_difference_geometric_max, u_difference_geometric_max],
                                            power = [k_difference_power_max, p_difference_power_max, u_difference_power_max],
                                            generalized = [k_difference_generalized_max, p_difference_generalized_max, u_difference_generalized_max],
                                            transmissibility = [k_difference_transmissibility_max, p_difference_transmissibility_max, u_difference_transmissibility_max])


fig = plt.figure()
#fig.suptitle("Field difference", fontsize=14)

ax = plt.subplot(3, 5, 1)
ax.set_title("Среднее степенное")
im_4 = ax.imshow(k_difference_generalized, extent = [0, 500, 0, 500])
cbar = colorbar(im_4)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.27, 0.5, "Поле абсолютной", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.20, 0.5, "проницаемости", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(3, 5, 2)
ax.set_title("Среднее гармоническое")
im_1 = ax.imshow(k_difference_harmonic, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)

ax = plt.subplot(3,5, 3)
ax.set_title("Среднее геометрическое")
im_2 = ax.imshow(k_difference_geometric, extent = [0, 500, 0, 500])
cbar = colorbar(im_2)
cbar.set_label(raw"$k, D$", fontsize = 14)

ax = plt.subplot(3, 5, 4)
ax.set_title("Усреднение мощности")
im_3 = ax.imshow(k_difference_power, extent = [0, 500, 0, 500])
cbar = colorbar(im_3)
cbar.set_label(raw"$k, D$", fontsize = 14)

ax = plt.subplot(3, 5, 5)
ax.set_title("Потоковый метод")
im_3 = ax.imshow(k_difference_transmissibility, extent = [0, 500, 0, 500])
cbar = colorbar(im_3)
cbar.set_label(raw"$k, D$", fontsize = 14)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-2, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(3,5,6)
im_7 = ax.imshow(p_difference_generalized, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_7)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
ax.text(-0.27, 0.5, "Поле", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.20, 0.5, "давления", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(3, 5, 7)
ax.set_title("Сред")
im_4 = ax.imshow(p_difference_harmonic, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_4)
cbar.set_label(raw"$P, MPa$", fontsize = 14)

ax = plt.subplot(3, 5, 8)
im_5 = ax.imshow(p_difference_geometric, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_5)
cbar.set_label(raw"$P, MPa$", fontsize = 14)

ax = plt.subplot(3, 5, 9)
im_6 = ax.imshow(p_difference_power, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_6)
cbar.set_label(raw"$P, MPa$", fontsize = 14)

ax = plt.subplot(3, 5, 10)
im_6 = ax.imshow(p_difference_transmissibility, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_6)
cbar.set_label(raw"$P, MPa$", fontsize = 14)

ax = plt.subplot(3, 5, 11)
im_10 = ax.imshow(u_difference_generalized, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_10)
cbar.set_label(raw"$U, m/s$", fontsize = 14)
ax.text(-0.27, 0.5, "Поле", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.20, 0.5, "скорости", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(3,5, 12)
im_7 = ax.imshow(u_difference_harmonic, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_7)
cbar.set_label(raw"$U, m/s$", fontsize = 14)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(3, 5, 13)
im_8 = ax.imshow(u_difference_geometric, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_8)
cbar.set_label(raw"$U, m/s$", fontsize = 14)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(3, 5, 14)
im_9 = ax.imshow(u_difference_power, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_9)
cbar.set_label(raw"$U, m/s$", fontsize = 14)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(3, 5, 15)
im_9 = ax.imshow(u_difference_transmissibility, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_9)
cbar.set_label(raw"$U, m/s$", fontsize = 14)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

plt.subplots_adjust(left=0.06, bottom=0.11, right=0.955, top=0.88, wspace=0.2, hspace=0.2)
plt.show()



c_0 = zeros(size(u_x, 1) - 19, size(u_y, 2) - 19)
for i = 1:size(u_x, 1) - 19
    c_0[i, 1] = 1
end
c_0_harmonic = zeros(size(u_x_harmonic, 1) - 3, size(u_y_harmonic, 2) - 3)
for i = 1:size(u_x_harmonic, 1) - 3
    c_0_harmonic[i, 1] = 1
end
c_0_geometric = zeros(size(u_x_geometric, 1) - 3, size(u_y_geometric, 2) - 3)
for i = 1:size(u_x_geometric, 1) - 3
    c_0_geometric[i, 1] = 1
end
c_0_power = zeros(size(u_x_power, 1) - 3, size(u_y_power, 2) - 3)
for i = 1:size(u_x_power, 1) - 3
    c_0_power[i, 1] = 1
end
c_0_generalized = zeros(size(u_x_generalized, 1) - 3, size(u_y_generalized, 2) - 3)
for i = 1:size(u_x_generalized, 1) - 3
    c_0_generalized[i, 1] = 1
end
c_0_transmissibility = zeros(size(u_x_transmissibility, 1) - 3, size(u_y_transmissibility, 2) - 3)
for i = 1:size(u_x_transmissibility, 1) - 3
    c_0_transmissibility[i, 1] = 1
end

coef = 0
t = 0
while coef < 0.05
   global t
   t += 1
   c = concentration_2d(c_0_harmonic, h_x*x, h_y*y, u_x_simple, u_y_simple, t)
   global coef = maximum(c[:, end])
end
coef_harmonic = 0
t_harmonic = 0
while coef_harmonic < 0
   global t_harmonic
   t_harmonic += 1
   c = concentration_2d(c_0_harmonic, h_x*x, h_y*y, u_x_harmonic, u_y_harmonic, t_harmonic)
   global coef_harmonic = maximum(c[:, end])
end
coef_geometric = 0
t_geometric = 0
while coef_geometric < 0.05
   global t_geometric
   t_geometric += 1
   c = concentration_2d(c_0_geometric, h_x*x, h_y*y, u_x_geometric, u_y_geometric, t_geometric)
   global coef_geometric = maximum(c[:, end])
end
coef_power = 0
t_power = 0
while coef_power < 0.05
   global t_power
   t_power += 1
   c = concentration_2d(c_0_power, h_x*x, h_y*y, u_x_power, u_y_power, t_power)
   global coef_power = maximum(c[:, end])
end
coef_generalized = 0
t_generalized = 0
while coef_generalized < 0.05
   global t_generalized
   t_generalized += 1
   c = concentration_2d(c_0_generalized, h_x*x, h_y*y, u_x_generalized, u_y_generalized, t_generalized)
   global coef_generalized = maximum(c[:, end])
end
coef_transmissibility = 0
t_transmissibility = 0
while coef_transmissibility < 0.05
   global t_transmissibility
   t_transmissibility += 1
   c = concentration_2d(c_0_transmissibility, h_x*x, h_y*y, u_x_transmissibility, u_y_transmissibility, t_transmissibility)
   global coef_transmissibility = maximum(c[:, end])
end
