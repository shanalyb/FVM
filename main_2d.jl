using Statistics, PyPlot, PyCall
@pyimport scipy.optimize as so
import PyPlot: quiver, colorbar, imshow
pygui(true)
include("gradient.jl")

N_x = 10; N_y = 10
L_x = 10; L_y = 10 #Расстояние (м)
P_left = 4; P_right = 2;
P_up = 4; P_down = 4;
h_x = L_x / N_x; h_y = L_y / N_y
dx = 10 / N_x; dy = 10 / N_y
#k = rand(Float64, N_x + 2, N_y + 2)
mu = 10^(-3)
h = 0.01
k = 0.5



a_old = zeros(N_x + 1, N_y + 1)
b_old = zeros(N_x + 1, N_y + 1)
A = zeros((N_x + 1) * (N_y + 1), (N_x + 1) * (N_y + 1))
b = zeros((N_x + 1) * (N_y + 1), 1)
p = zeros(N_x + 1, N_y + 1)
p_min = zeros(N_y + 1, 1)
p_max = zeros(N_y + 1, 1)
p_ex = zeros(1, (N_x + 1) * (N_y + 1))

for i = 1:N_x+1, j = 1:N_x+1
    a_old[i, j] = (k / h_x)
end

for i = 1:N_y+1, j = 1:N_y+1
    b_old[i, j] = (k / h_y)
end

#for i = 2:N_x+2, j = 1:N_x+1
#    a_old[i-1, j] = (h_y / (k[i, j] - k[i-1, j])) * (log(k[i, j] / k[i-1, j]))
#end

#for i = 1:N_y+1, j = 2:N_y+2
#    b_old[i, j-1] = (h_x / (k[i, j] - k[i, j-1])) * (log(k[i, j] / k[i, j-1]))
#end1


for m = 1:(N_y-1)
    A[m, m+1] = 1
    A[N_y-1+m, N_x*(N_y+1)+m+1] = 1
    b[m, 1] = P_left
    b[N_y-1+m, 1] = P_right
end

for m = 0:N_x
    A[2*(N_y-1)+(m+1), m*(N_y+1)+1] = 1
    A[2*(N_y-1)+(N_x+1)+(m+1), m*(N_y+1)+(N_y+1)] = 1
    #b[2*(N_y-1)+(m+1), 1] = P_up
    #b[2*(N_y-1)+(N_x+1)+(m+1), 1] = P_down
    A[2*(N_y-1)+(m+1), m*(N_y+1)+1+1] = -1
    A[2*(N_y-1)+(N_x+1)+(m+1), m*(N_y+1)+(N_y+1)-1] = -1
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

fig = imshow(p, extent=[0,L_x,0,L_y], interpolation = "gaussian")
#fig.set_xlabel("X")
cbar = colorbar(fig)
cbar.set_label("P, MPa", fontsize=14)

u_x = - (k / mu) * grad2D_x(p, h_x)
u_y = - (k / mu) * grad2D_y(p, h_y)

#x = zeros(L_x, L_y), y = zeros(L_x, L_y)

#function numpy_mgrid(dim1,dim2)
#    X = [i for i in 1:dim1, j in 1:dim2]
#    Y = [j for i in 1:dim1, j in 1:dim2]
#    return X,Y
#end

#x, y = numpy_mgrid(L_x, L_y)
#x, y = 0.5*x, 0.5*y
#x = linspace(0,1,L_x)
#y = linspace(1,0,L_y)
#x = range(0, stop=1, length=L_x)
#y = range(1, stop=0, length=L_y)
quiver(u_x, u_y)
