using PyCall, Statistics
@pyimport matplotlib.animation as anim
using PyPlot
pygui(true)
include("gradient.jl")
include("concentration_2d.jl")

N_x = 10;
N_y = 10;
L_x = 10; #Расстояние (м)
L_y = 10; #Расстояние (м)
P_left = 4;
P_right = 2;
P_up = 4;
P_down = 4;
h_x = L_x / N_x;
h_y = L_y / N_y;
dx = 10 / N_x;
dy = 10 / N_y;
k = rand(Float64, N_x + 2, N_y + 2)
mu = 10^(-3)
h = 0.01
#k = 0.5



a_old = zeros(N_x + 1, N_y + 1)
b_old = zeros(N_x + 1, N_y + 1)
A = zeros((N_x + 1) * (N_y + 1), (N_x + 1) * (N_y + 1))
b = zeros((N_x + 1) * (N_y + 1), 1)
p = zeros(N_x + 1, N_y + 1)
p_min = zeros(N_y + 1, 1)
p_max = zeros(N_y + 1, 1)
p_ex = zeros(1, (N_x + 1) * (N_y + 1))

#for i = 1:N_x+1, j = 1:N_x+1
#    a_old[i, j] = (k / h_x)
#end

#for i = 1:N_y+1, j = 1:N_y+1
#    b_old[i, j] = (k / h_y)
#end

#Неоднородная среда_____________
for i = 2:N_x+2, j = 1:N_x+1
    a_old[i-1, j] = (h_y / (k[i, j] - k[i-1, j])) * (log(k[i, j] / k[i-1, j]))
end

for i = 1:N_y+1, j = 2:N_y+2
    b_old[i, j-1] = (h_x / (k[i, j] - k[i, j-1])) * (log(k[i, j] / k[i, j-1]))
end
#______________________________

for m = 1:(N_y-1)
    A[m, m+1] = 1
    A[N_y-1+m, N_x * (N_y + 1) + m + 1] = 1
    b[m, 1] = P_left
    b[N_y-1+m, 1] = P_right
end

for m = 0:N_x
    A[2*(N_y-1)+(m+1), m*(N_y+1)+1] = 1
    A[2 * (N_y - 1) + (N_x + 1) + (m + 1), m*(N_y+1)+(N_y+1)] = 1
    #b[2*(N_y-1)+(m+1), 1] = P_up
    #b[2*(N_y-1)+(N_x+1)+(m+1), 1] = P_down
    A[2*(N_y-1)+(m+1), m * (N_y + 1) + 1 + 1] = -1
    A[2 * (N_y - 1) + (N_x + 1) + (m + 1), m*(N_y+1)+(N_y+1)-1] = -1
end

m = 2 * (N_x + N_y)
for i = 1:N_x-1, j = 1:N_y-1
    global m = m + 1
    A[m, (i-1)*(N_y+1)+(j+1)] = -a_old[i, j]
    A[m, i*(N_y+1)+(j+1)] = a_old[i, j] + a_old[i+1, j] + b_old[i, j] +
                            b_old[i, j+1]
    A[m, (i+1)*(N_y+1)+(j+1)] = -a_old[i+1, j]
    A[m, i*(N_y+1)+j] = -b_old[i, j]
    A[m, i*(N_y+1)+(j+2)] = -b_old[i, j+1]
end

A_new = inv(A)

p_old = A_new * b

n = 0
for j = 1:N_y+1, i = 1:N_x+1
    global n = n + 1
    p[i, j] = p_old[n]
end

#fig = imshow(p, extent = [0, L_x, 0, L_y], interpolation = "gaussian")
#cbar = colorbar(fig)
#cbar.set_label("P, MPa", fontsize = 14)

#Неоднородная среда_____________
k_x = zeros(N_x, N_y + 1)
for i = 1:N_x, j = 1:N_y+1
    k_x[i, j] = k[i, j]
end
k_y = zeros(N_x + 1, N_y)
for i = 1:N_x+1, j = 1:N_y
    k_y[i, j] = k[i, j]
end
#_______________________________
u_x = zeros(N_x + 1, N_y)
u_y = zeros(N_x, N_y + 1)
for i = 1:N_x+1, j = 1:N_y
    u_x[i, j] = -(k[i, j] / mu) * grad2D_x(p, h_x)[i, j] * 10^(-6)
end
for i = 1:N_x, j = 1:N_y+1
    u_y[i, j] = -(k[i, j] / mu) * grad2D_y(p, h_y)[i, j] * 10^(-6)
end

#quiver(u_x, u_y)

#Концентрация
c_0 = zeros(N_x + 1, N_y + 1)
for i = 1:N_x+1
    c_0[i, 1] = 1
end
#Неоднородная среда_____________
#c_0 = zeros(N_x, N_y)
#for i = 1:N_x
#    c_0[i, 1] = 1
#end
#u_x_si_inhomogenes = zeros(N_x, N_y-1)
#for i = 1:N_x, j = N_y-1
#    u_x_si_inhomogenes[i, j] = u_x_si[i, j]
#end
#u_y_si_inhomogenes = zeros(N_x-1, N_y)
#for i = 1:N_x-1, j = 1:N_y
#    u_y_si_inhomogenes[i, j] = u_y_si[i, j]
#end
#_______________________________

t = 200
c = concentration_2d(c_0, N_x, N_y, h_x, h_y, u_x, u_y, t)

fig = imshow(c, extent = [0, L_x + 1, 0, L_y + 1], interpolation = "gaussian")
cbar = colorbar(fig)
cbar.set_label("Concentration", fontsize = 14)
