using PyCall
import PyPlot: quiver
using PyPlot
pygui(true)
include("main_2d.jl")

x = range(0, stop=480, length=49)
y = range(480, stop=0, length=49)

x_upscaling = range(0, stop=500, length=16)
y_upscaling = range(500, stop=0, length=16)

fig = plt.figure()
#fig.suptitle("Сравнение методов апскейлинга", fontsize=14)

#Проницаемости
ax = plt.subplot(6,3,1)
ax.set_title("Поле проницаемости", y=1.1)
im_1 = ax.imshow(k, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Мелкая", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "сетка", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(6,3,4)
im_1 = ax.imshow(k_generalized, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Среднее", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "степенное", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(6,3,7)
im_1 = ax.imshow(k_harmonic, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Среднее", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "гармоническое", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(6,3,10)
im_1 = ax.imshow(k_geometric, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Среднее", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "геометрическое", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(6,3,13)
im_1 = ax.imshow(k_power, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Усреднение", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "мощности", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

ax = plt.subplot(6,3,16)
im_1 = ax.imshow(k_transmissibility, extent = [0, 500, 0, 500])
cbar = colorbar(im_1)
cbar.set_label(raw"$k, D$", fontsize = 14)
ax.text(-0.47, 0.5, "Потоковый", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)
ax.text(-0.35, 0.5, "метод", horizontalalignment="right", verticalalignment="center", rotation="vertical", transform=ax.transAxes, size = 12)

#_____________________

#Давления
ax = plt.subplot(6,3,2)
ax.set_title("Поле давления и скорости", y=1.1)
im_2 = ax.imshow(p, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_2)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x, y, u_x, u_y)

ax = plt.subplot(6,3,5)
im_3 =ax.imshow(p_generalized, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x_upscaling, y_upscaling, u_x_harmonic, u_y_harmonic)

ax = plt.subplot(6,3,8)
im_3 =ax.imshow(p_harmonic, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x_upscaling, y_upscaling, u_x_harmonic, u_y_harmonic)

ax = plt.subplot(6,3,11)
im_3 =ax.imshow(p_geometric, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x_upscaling, y_upscaling, u_x_harmonic, u_y_harmonic)

ax = plt.subplot(6,3,14)
im_3 =ax.imshow(p_power, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x_upscaling, y_upscaling, u_x_harmonic, u_y_harmonic)

ax = plt.subplot(6,3,17)
im_3 =ax.imshow(p_transmissibility, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_3)
cbar.set_label(raw"$P, MPa$", fontsize = 14)
quiver(x_upscaling, y_upscaling, u_x_harmonic, u_y_harmonic)
#_____________________

#Отношение скоростей
ax = plt.subplot(6,3,3)
ax.set_title("Отношение u_y/u_x", y=1.1)
im_4 =ax.imshow(koef, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(6,3,6)
im_4 =ax.imshow(koef_generalized, extent = [0, 500, 0, 500], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(6,3,9)
im_4 =ax.imshow(koef_harmonic, extent = [0, size(koef, 1)*10, 0, size(koef, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(6,3,12)
im_4 =ax.imshow(koef_geometric, extent = [0, size(koef, 1)*10, 0, size(koef, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(6,3,15)
im_4 =ax.imshow(koef_power, extent = [0, size(koef, 1)*10, 0, size(koef, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)

ax = plt.subplot(6,3,18)
im_4 =ax.imshow(koef_transmissibility, extent = [0, size(koef, 1)*10, 0, size(koef, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_4)
formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
cbar.ax.yaxis.set_major_formatter(formatter)
#_____________________


plt.subplots_adjust(left=0.265, bottom=0.025, right=0.73, top=0.955, wspace=0.0, hspace=0.225)
plt.show()


fig = plt.figure()
fig.suptitle("Концентрация примеси", fontsize=14)
#Концентрация
ax = plt.subplot(1,6,1)
ax.set_title("Мелкая сетка", y=1.1)
im_3 =ax.imshow(c, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)

ax = plt.subplot(1,6,2)
ax.set_title("Среднее степенное", y=1.1)
im_3 =ax.imshow(c_generalized, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)

ax = plt.subplot(1,6,3)
ax.set_title("Среднее гармоническое", y=1.1)
im_3 =ax.imshow(c_harmonic, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)

ax = plt.subplot(1,6,4)
ax.set_title("Среднее геометрическое", y=1.1)
im_3 =ax.imshow(c_geometric, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)

ax = plt.subplot(1,6,5)
ax.set_title("Усреднение мощности", y=1.1)
im_3 =ax.imshow(c_power, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)

ax = plt.subplot(1,6,6)
ax.set_title("Потоковый метод", y=1.1)
im_3 =ax.imshow(c_transmissibility, extent = [0, size(c, 1)*10, 0, size(c, 2)*10], interpolation = "gaussian")
cbar = colorbar(im_3)
#_____________________
plt.subplots_adjust(left=0.0, bottom=0.65, right=0.995, top=0.87, wspace=0.0, hspace=0.3)
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

#ax = plt.subplot(1,6,4)
#ax.set_title("generalized average")
#im_1 = ax.imshow(k_generalized, extent = [0, size(k_geometric, 1), 0, size(k_geometric, 2)])


#plt.show()
