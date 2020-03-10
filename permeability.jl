using PyCall, Statistics, PyPlot, LinearAlgebra
pygui(true)
N_x = 10
N_y = 10
r_max = 3
sigma = 1
k_x = 1
k_y = 1

k = zeros(N_y + 2, N_x + 2)
k_rand = zeros(N_y + 2, N_x + 2)
R = zeros((N_x+2)*(N_y+2), (N_x+2)*(N_y+2))
R_new = zeros((N_x+2)*(N_y+2), (N_x+2)*(N_y+2))
for i = 1:(N_x+2)*(N_y+2)
    for j = 1:(N_x+2)*(N_y+2)
        x = abs(i%(N_x+2)-j%(N_y+2))/k_x
        y = abs(div(i, N_x+2)-div(j, N_y+2))/k_y
        r = sqrt(x^2+y^2)
        if r <= r_max
            R[j, i] = ((sigma^2) / 12)*((r_max-r)/r_max)
        end
    end
end
fig = imshow(R, extent = [0, (N_x+1)*(N_y+1), 0, (N_x+1)*(N_y+1)], interpolation = "gaussian")
cbar = colorbar(fig)
R_new = sqrt(R)
