using PyCall, DataFrames
import PyPlot: quiver
using PyPlot
pygui(true)
include("main_2d.jl")
include("upscaling.jl")
k_simple = simple_average(k, x, y)
p_simple = simple_average(p, x, y)
u_simple = simple_average(u_summ, x, y)


function equalize_the_sizes(p, p_small)
    p_new = zeros(size(p_small, 1), size(p_small, 2))
    for i = 1:size(p_small, 1), j = 1:size(p_small, 2)
        p_new[i, j] = p[i, j]
    end
    return p_new
end

p_harmonic = equalize_the_sizes(p_harmonic, p_simple)
p_geometric = equalize_the_sizes(p_geometric, p_simple)
p_power = equalize_the_sizes(p_power, p_simple)
p_generalized = equalize_the_sizes(p_generalized, p_simple)

u_simple = equalize_the_sizes(u_simple, u_geometric)

k_difference_harmonic = abs.(k_simple - k_harmonic)
k_difference_geometric = abs.(k_simple - k_geometric)
k_difference_power = abs.(k_simple - k_power)
k_difference_generalized = abs.(k_simple - k_generalized)

p_difference_harmonic = abs.(p_simple - p_harmonic)
p_difference_geometric = abs.(p_simple - p_geometric)
p_difference_power = abs.(p_simple - p_power)
p_difference_generalized = abs.(p_simple - p_generalized)

u_difference_harmonic = abs.(u_simple - u_harmonic)
u_difference_geometric = abs.(u_simple - u_geometric)
u_difference_power = abs.(u_simple - u_power)
u_difference_generalized = abs.(u_simple - u_generalized)

#Сравнение отклонений
k_difference_harmonic_max = maximum(k_difference_harmonic)
k_difference_geometric_max = maximum(k_difference_geometric)
k_difference_power_max = maximum(k_difference_power)
k_difference_generalized_max = maximum(k_difference_generalized)

p_difference_harmonic_max = maximum(p_difference_harmonic)
p_difference_geometric_max = maximum(p_difference_geometric)
p_difference_power_max = maximum(p_difference_power)
p_difference_generalized_max = maximum(p_difference_generalized)

u_difference_harmonic_max = maximum(u_difference_harmonic)
u_difference_geometric_max = maximum(u_difference_geometric)
u_difference_power_max = maximum(u_difference_power)
u_difference_generalized_max = maximum(u_difference_generalized)

df = DataFrame(parameters = ["k", "p", "u"], harmonic = [k_difference_harmonic_max, p_difference_harmonic_max, u_difference_harmonic_max], geometric = [k_difference_geometric_max, p_difference_geometric_max, u_difference_geometric_max], power = [k_difference_power_max, p_difference_power_max, u_difference_power_max], generalized = [k_difference_generalized_max, p_difference_generalized_max, u_difference_generalized_max])


fig = plt.figure()
fig.suptitle("Field difference", fontsize=14)

ax = plt.subplot(3, 4, 1)
ax.set_title("difference of concentration fields (harmonic)")
im_1 = ax.imshow(k_difference_harmonic, extent = [0, size(k_difference_harmonic, 1), 0, size(k_difference_harmonic, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(3,4, 2)
ax.set_title("difference of concentration fields (geometric)")
im_2 = ax.imshow(k_difference_geometric, extent = [0, size(k_difference_geometric, 1), 0, size(k_difference_geometric, 2)])
cbar = colorbar(im_2)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(3, 4, 3)
ax.set_title("difference of concentration fields (power)")
im_3 = ax.imshow(k_difference_power, extent = [0, size(k_difference_power, 1), 0, size(k_difference_power, 2)])
cbar = colorbar(im_3)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(3, 4, 4)
ax.set_title("difference of concentration fields (generalized)")
im_4 = ax.imshow(k_difference_generalized, extent = [0, size(k_difference_generalized, 1), 0, size(k_difference_generalized, 2)])
cbar = colorbar(im_4)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(3, 4, 5)
ax.set_title("difference of pressure fields (harmonic)")
im_4 = ax.imshow(p_difference_harmonic, extent = [0, size(p_difference_harmonic, 1)-1, 0, size(p_difference_harmonic, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_4)
cbar.set_label("P, MPa", fontsize = 14)

ax = plt.subplot(3, 4, 6)
ax.set_title("difference of pressure fields (geometric)")
im_5 = ax.imshow(p_difference_geometric, extent = [0, size(p_difference_geometric, 1)-1, 0, size(p_difference_geometric, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_5)
cbar.set_label("P, MPa", fontsize = 14)

ax = plt.subplot(3, 4, 7)
ax.set_title("difference of pressure fields (power)")
im_6 = ax.imshow(p_difference_power, extent = [0, size(p_difference_power, 1)-1, 0, size(p_difference_power, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_6)
cbar.set_label("P, MPa", fontsize = 14)

ax = plt.subplot(3,4, 8)
ax.set_title("difference of pressure fields (generalized)")
im_7 = ax.imshow(p_difference_generalized, extent = [0, size(p_difference_generalized, 1)-1, 0, size(p_difference_generalized, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_7)
cbar.set_label("P, MPa", fontsize = 14)

ax = plt.subplot(3,4, 9)
ax.set_title("difference of velocity fields (harmonic)")
im_7 = ax.imshow(u_difference_harmonic, extent = [0, size(u_difference_harmonic, 1), 0, size(u_difference_harmonic, 2)], interpolation = "gaussian")
cbar = colorbar(im_7)
cbar.set_label("U, m/s", fontsize = 14)

ax = plt.subplot(3, 4, 10)
ax.set_title("difference of velocity fields (geometric)")
im_8 = ax.imshow(u_difference_geometric, extent = [0, size(u_difference_geometric, 1), 0, size(u_difference_geometric, 2)], interpolation = "gaussian")
cbar = colorbar(im_8)
cbar.set_label("U, m/s", fontsize = 14)

ax = plt.subplot(3, 4, 11)
ax.set_title("difference of velocity fields (power)")
im_9 = ax.imshow(u_difference_power, extent = [0, size(u_difference_power, 1), 0, size(u_difference_power, 2)], interpolation = "gaussian")
cbar = colorbar(im_9)
cbar.set_label("U, m/s", fontsize = 14)

ax = plt.subplot(3, 4, 12)
ax.set_title("difference of velocity fields (generalized)")
im_10 = ax.imshow(u_difference_generalized, extent = [0, size(u_difference_generalized, 1), 0, size(u_difference_generalized, 2)], interpolation = "gaussian")
cbar = colorbar(im_10)
cbar.set_label("U, m/s", fontsize = 14)

plt.show()
