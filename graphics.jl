using PyCall
import PyPlot: quiver
using PyPlot
pygui(true)
include("main_2d.jl")

fig = plt.figure()
fig.suptitle("Original and power average", fontsize=14)

ax = plt.subplot("241")
ax.set_title("concentration field")
im_1 = ax.imshow(k, extent = [0, size(k, 1), 0, size(k, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot("242")
ax.set_title("pressure and velocity")
im_2 = ax.imshow(p, extent = [0, size(p, 1)-1, 0, size(p, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_2)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x, u_y)

ax = plt.subplot("243")
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c, extent = [0, size(c, 1), 0, size(c, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot("244")
ax.set_title("u_y/u_x")
im_4 =ax.imshow(koef, extent = [0, size(koef, 1), 0, size(koef, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)

ax = plt.subplot("245")
ax.set_title("concentration field")
im_1 = ax.imshow(k_geometric, extent = [0, size(k_harmonic, 1), 0, size(k_harmonic, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot("246")
ax.set_title("pressure and velocity")
im_3 =ax.imshow(p_geometric, extent = [0, size(p_harmonic, 1)-1, 0, size(p_harmonic, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x_harmonic, u_y_harmonic)

ax = plt.subplot("247")
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c_geometric, extent = [0, size(c_harmonic, 1), 0, size(c_harmonic, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot("248")
ax.set_title("u_y/u_x")
im_4 =ax.imshow(koef_geometric, extent = [0, size(koef_harmonic, 1), 0, size(koef_harmonic, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)

plt.show()
