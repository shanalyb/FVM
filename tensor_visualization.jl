using Makie, Images, ImageFiltering, ImageView, ImageSegmentation, LinearAlgebra

function draw_tensor(K, filename)
    # Utility to split an array-of-length-3-arrays into x,y,z components
    split_coords(ps) = getindex.(ps,1), getindex.(ps,2), getindex.(ps,3)

    # Utility to convert vectors to colors
    to_rgb(v) = RGBf0(v[1], v[2], v[3])

    # Convert spherical coordinates to Cartesian coordinates.
    # Note that θ,ϕ are named according to the "physicist convention"
    # Note that using Vec3f0 here is efficient, but you could just use a
    # normal array and everything would still work.
    spherical_to_cartesian(r,θ,ϕ) = Vec3f0(r*sin(θ)*cos(ϕ), r*sin(θ)*sin(ϕ), r*cos(θ))

    θ = range(0,pi,length=20)
    ϕ = permutedims(range(0,2pi,length=30))

    # A parametric unit sphere
    v = spherical_to_cartesian.(1.0, θ, ϕ)

    # One way of "seeing how K acts in each direction" is by taking the norm of
    # K * v and using this as the radius in terms of the original angular
    # coordinates.
    #r = norm.(Ref(K) .* v)
    #ps = spherical_to_cartesian.(r, θ, ϕ)

    # Another way is by seeing how it distorts the unit sphere by computing K*v for
    # every point v on the sphere:
    ps = Ref(K) .* v

    # Ultimately both of the above techniques only present part of the information
    # about K in the position because K contains more degrees of freedom than we
    # can represent on an ellipse.

    x,y,z = split_coords(ps)
    scene = surface(x,y,z, color=to_rgb.(v))

    x,y,z = split_coords(1.003.*ps) # Scaling hack to avoid some z-fighting
    wireframe!(scene, x,y,z)

    Makie.save(filename, scene)
end

 #K =[ 0.070087    -0.0386396   0.00209088;
#     -0.0441344    0.105735    0.00012994;
#     -0.00553273  -0.00292807  0.0896236]

#draw_tensor(K, "example.png")
