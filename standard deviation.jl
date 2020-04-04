using Statistics
include("upscaling.jl")
include("main_2d.jl")

#Отклонение
#--------------------------------
k_otklon = std(k)
k_otklon_harmonic = std(k_harmonic)
k_otklon_geometric = std(k_geometric)
k_otklon_power = std(k_power)

p_otklon = std(p)
p_otklon_harmonic = std(p_harmonic)
p_otklon_geometric = std(p_geometric)
p_otklon_power = std(p_power)

u_otklon = std(u_summ)
u_otklon_harmonic = std(u_harmonic)
u_otklon_geometric = std(u_geometric)
u_otklon_power = std(u_power)

#println("ср.квадратичн. отклон. = ", k_otklon)
#println("ср.квадратичн. отклон. harmonic = ", k_otklon_harmonic)
#println("ср.квадратичн. отклон. geometric = ", k_otklon_geometric)
#println("ср.квадратичн. отклон. power = ", k_otklon_power)
#--------------------------------

fig = plt.figure()
color_rectangle = ["b", "g", "r", "c"]

ax = plt.subplot("131")
ax.set_title("standard deviation of the concentration")
x = [1, 2, 3, 4]
y = [k_otklon, k_otklon_harmonic, k_otklon_geometric, k_otklon_power]
ax.bar(x, y, color = color_rectangle, width = 0.8)
plt.xticks(x, ["original","harmonic","geometric","power"])
ax.set_facecolor("seashell")
ax.set_ylabel("k, D", fontsize = 14)

ax = plt.subplot("132")
ax.set_title("the standard deviation of the pressure")
x = [1, 2, 3, 4]
y = [p_otklon, p_otklon_harmonic, p_otklon_geometric, p_otklon_power]
ax.bar(x, y, color = color_rectangle, width = 0.8)
plt.xticks(x, ["original","harmonic","geometric","power"])
ax.set_facecolor("seashell")
ax.set_ylabel("P, MPa", fontsize = 14)

ax = plt.subplot("133")
ax.set_title("standard deviation of the velocity")
x = [1, 2, 3, 4]
y = [u_otklon, u_otklon_harmonic, u_otklon_geometric, u_otklon_power]
ax.bar(x, y, color = color_rectangle, width = 0.8)
plt.xticks(x, ["original","harmonic","geometric","power"])
ax.set_facecolor("seashell")
ax.set_ylabel("U, m/s", fontsize = 14)

formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3, 2))
ax.yaxis.set_major_formatter(formatter)

fig.set_facecolor("floralwhite")

plt.show()
