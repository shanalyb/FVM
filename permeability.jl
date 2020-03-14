using PyCall, Statistics, PyPlot, LinearAlgebra
pygui(true)
N_x = 7
N_y = 7
r_max = 1
sigma = 1
k_x = 1
k_y = 1
a = rand(Float64, N_x+1, N_y+1)
a_new = reshape(a, ((N_x+1)*(N_y+1), 1))

h = sigma^2/12
CovMatrix = zeros((N_x+1)*(N_y+1), (N_x+1)*(N_y+1))
for i = 0:N_x, j = 0:N_y, l = 0:N_x, m = 0:N_y
    k = i*(N_y+1)+(j+1)
    n = l*(N_y+1)+(m+1)
    r = (((i-l)/k_x)^2+((j-m)/k_y)^2)^0.5
    if r <= r_max
            CovMatrix[k, n] = h*(r_max-r)/r_max
        else
            CovMatrix[k, n] = 0
    end
end
R = sqrt(CovMatrix)
subplot(1, 2, 1)
fig = imshow(R, extent = [0, (N_x+1)*(N_y+1), 0, (N_x+1)*(N_y+1)], interpolation = "gaussian")
cbar = colorbar(fig)

nu = R*a_new
nu_matrix = reshape(nu, (N_x+1, N_y+1))
subplot(1, 2, 2)
fig = imshow(nu_matrix, extent = [0, (N_x+1)*(N_y+1), 0, (N_x+1)*(N_y+1)])
cbar = colorbar(fig)
#R[j, i] = ((sigma^2) / 12)*((r_max-r)/r_max)
