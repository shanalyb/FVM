using PyCall
import PyPlot: quiver
using PyPlot
pygui(true)
include("main_2d.jl")

fig = plt.figure()
fig.suptitle("Original and power average", fontsize=14)

ax = plt.subplot(4,4,1)
ax.set_title("the permeability field")
im_1 = ax.imshow(k, extent = [0, size(k, 1), 0, size(k, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(4,4,2)
ax.set_title("pressure and velocity")
im_2 = ax.imshow(p, extent = [0, size(p, 1)-1, 0, size(p, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_2)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x, u_y)

ax = plt.subplot(4,4,3)
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c, extent = [0, size(c, 1), 0, size(c, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot(4,4,4)
ax.set_title("u_y/u_x")
im_4 =ax.imshow(u_summ, extent = [0, size(koef, 1), 0, size(koef, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)


ax = plt.subplot(4,4,5)
ax.set_title("the permeability field")
im_1 = ax.imshow(k_generalized, extent = [0, size(k_generalized, 1), 0, size(k_generalized, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(4,4,6)
ax.set_title("pressure and velocity")
im_3 =ax.imshow(p_generalized, extent = [0, size(p_generalized, 1)-1, 0, size(p_generalized, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x_harmonic, u_y_harmonic)

ax = plt.subplot(4,4,7)
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c_generalized, extent = [0, size(c_generalized, 1), 0, size(c_generalized, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot(4,4,8)
ax.set_title("u_y/u_x")
im_4 =ax.imshow(koef_generalized, extent = [0, size(koef_generalized, 1), 0, size(koef_generalized, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)



ax = plt.subplot(4,4,9)
ax.set_title("the permeability field")
im_1 = ax.imshow(k_harmonic, extent = [0, size(k_harmonic, 1), 0, size(k_harmonic, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(4,4,10)
ax.set_title("pressure and velocity")
im_3 =ax.imshow(p_harmonic, extent = [0, size(p_harmonic, 1)-1, 0, size(p_harmonic, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x_harmonic, u_y_harmonic)

ax = plt.subplot(4,4,11)
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c_harmonic, extent = [0, size(c_harmonic, 1), 0, size(c_harmonic, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot(4,4,12)
ax.set_title("u_y/u_x")
im_4 =ax.imshow(koef_harmonic, extent = [0, size(koef_harmonic, 1), 0, size(koef_harmonic, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)



ax = plt.subplot(4,4,13)
ax.set_title("the permeability field")
im_1 = ax.imshow(k_geometric, extent = [0, size(k_geometric, 1), 0, size(k_geometric, 2)])
cbar = colorbar(im_1)
cbar.set_label("k, D", fontsize = 14)

ax = plt.subplot(4,4,14)
ax.set_title("pressure and velocity")
im_3 =ax.imshow(p_geometric, extent = [0, size(p_geometric, 1)-1, 0, size(p_geometric, 2)-1], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("P, MPa", fontsize = 14)
quiver(u_x_harmonic, u_y_harmonic)

ax = plt.subplot(4,4,15)
ax.set_title("concentration of impurity")
im_3 =ax.imshow(c_geometric, extent = [0, size(c_geometric, 1), 0, size(c_geometric, 2)], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label("Concentration", fontsize = 14)

ax = plt.subplot(4,4,16)
ax.set_title("u_y/u_x")
im_4 =ax.imshow(koef_geometric, extent = [0, size(koef_geometric, 1), 0, size(koef_geometric, 2)], interpolation = "gaussian")
cbar = colorbar(im_4)

plt.subplots_adjust(wspace=0, hspace=0.3)
plt.show()


#fig = plt.figure()
#fig.suptitle("Original and power average", fontsize=14)

#ax = plt.subplot(1,4,1)
#ax.set_title("original")
#im_1 = ax.imshow(k, extent = [0, size(k, 1), 0, size(k, 2)])

#ax = plt.subplot(1,4,2)
#ax.set_title("harmonic average")
#im_1 = ax.imshow(k_harmonic, extent = [0, size(k_generalized, 1), 0, size(k_generalized, 2)])

#ax = plt.subplot(1,4,3)
#ax.set_title("geometric average")
#im_1 = ax.imshow(k_geometric, extent = [0, size(k_harmonic, 1), 0, size(k_harmonic, 2)])

#ax = plt.subplot(1,4,4)
#ax.set_title("generalized average")
#im_1 = ax.imshow(k_generalized, extent = [0, size(k_geometric, 1), 0, size(k_geometric, 2)])


#plt.show()
