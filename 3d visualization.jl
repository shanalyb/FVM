using Makie
using LinearAlgebra
using StaticArrays

# Utility to split an array-of-length-3-arrays into x,y,z components
split_coords(ps) = getindex.(ps,1), getindex.(ps,2), getindex.(ps,3)

# Spherical Fibonacci numbers for Quasi-random distribution on the sphere.
# See "Spherical Fibonacci Point Sets for Illumination Integrals" by Marques et al.,
# doi:10.1111/cgf.12190
function spherical_fib(j, N)
θ = acos(1 - 2*j/N)
ϕ = 2*j*π/MathConstants.φ
(θ,ϕ)
end

# Convert spherical coordinates to Cartesian coordinates.
# Note that θ,ϕ are named according to the "physicist convention"
# Note that using Vec3f0 here is efficient, but you could just use a
# normal array and everything would still work.
spherical_to_cartesian(r,θ,ϕ) = Vec3f0(r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ))

scene = Scene()
for x in -3:3:3, y = -3:3:3
origin = Point3f0(x,y,0)
K = @SArray randn(3,3)
# Plot a parametric unit sphere
θ = range(0,pi,length=100)
ϕ = permutedims(range(0,2pi,length=100))
ps = spherical_to_cartesian.(1.0, θ, ϕ)
surface!(scene, split_coords(Ref(origin) .+ ps)..., color=norm.(Ref(K) .* ps))

N = 500
vs = map(j->spherical_to_cartesian(1, spherical_fib(j,N)...), 1:N)
Kvs = Ref(K) .* vs
arrows!(scene, Ref(origin) .+ vs, 0.1 .* Kvs, arrowsize=0.05, linewidth=3, lengthscale=3, arrowcolor=norm.(Kvs))
end
scene

# Makie.save("tensor.png", scene)
