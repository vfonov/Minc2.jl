# for Affine transforms
using LinearAlgebra

# for grid transforms
using Interpolations

"""
Affine transform
"""
mutable struct AffineTransform
    mat::Matrix{Float64}
end

function AffineTransform()
    AffineTransform([1.0 0.0 0.0 0.0;0.0 1.0 0.0 0.0;0.0 0.0 1.0 0.0;0.0 0.0 0.0 1.0])
end


"""
Dense vector field transform (grid transform)
"""
mutable struct GridTransform
    voxel_to_world::AffineTransform
    world_to_voxel::AffineTransform
    inverse::Bool
    vector_field::Array{Float64, 4}
    itp_vector_field
    function GridTransform(voxel_to_world::AffineTransform,
        inverse::Bool,
        vector_field::Array{Float64, 4})
        new(voxel_to_world,inv(voxel_to_world),inverse,vector_field,
            extrapolate(interpolate(vector_field, (NoInterp(),BSpline(Linear()),BSpline(Linear()),BSpline(Linear()))),Flat()))
    end
end

function GridTransform()
    GridTransform(
        AffineTransform(),
        false,
        zeros(3,3,3,3)
    )
end

"""
AnyTransform
"""
AnyTransform=Union{AffineTransform,GridTransform}


"""
Invert AffineTransform transform
"""
function inv(t::AffineTransform)::AffineTransform
    AffineTransform(Base.inv(t.mat))
end

"""
Invert GridTransform transform
"""
function inv(t::GridTransform)::GridTransform
    GridTransform(t.voxel_to_world, !t.inverse, t.vector_field)
end

"""
Invert concatenated transform
"""
function inv(t::Vector{AnyTransform})::Vector{AnyTransform}
    [inv(i) for i in reverse(t)]
end

"""
Apply affine transform
"""
function transform_point(tfm::AffineTransform,p::Vector{Float64})::Vector{Float64}
    (p' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
end

"""
Apply grid transform
# TODO: handle inversion flag?
"""
function transform_point(tfm::GridTransform, p::Vector{Float64})::Vector{Float64}
    #(p' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
    # convert to voxel coords, add 1 to get index
    v = transform_point(tfm.world_to_voxel, p) .+ 1.0
    # interpolate displacement vector
    d=[tfm.itp_vector_field(1,v...),
       tfm.itp_vector_field(2,v...),
       tfm.itp_vector_field(3,v...) ]
    # resulting deformation
    
    return p+d
end


"""
Apply concatenated transform
"""
function transform_point(tfm::Vector{AnyTransform}, p::Vector{Float64})::Vector{Float64}
    for t in tfm
        p=transform_point(t,p)
    end
    return p
end


"""
Apply affine transform to CartesianIndices
"""
function transform_point(tfm::AffineTransform,p::CartesianIndex{3})::Vector{Float64}
    ( [p[1]-1.0,p[2]-1.0,p[3]-1.0]' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
end


# helper 
Base.show(io::IO, z::GridTransform) = print(io, "GridTransform:", size(z.vector_field))
