using LinearAlgebra

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
    world_to_voxel::AffineTransform
    voxel_to_world::AffineTransform
    inverse::Bool
    vector_field::Array{Float64, 4}
end

function GridTransform()
    GridTransform(
        AffineTransform(),
        AffineTransform(),
        false,
        zeros(1,1,1,3)
    )
end

"""
AnyTransform
"""
AnyTransform=Union{AffineTransform,GridTransform}


"""
Invert transform
"""
function inv(t::AffineTransform)::AffineTransform
    AffineTransform(Base.inv(t.mat))
end

"""
Apply affine transform
"""
function transform_point(tfm::AffineTransform,p::Vector{Float64})::Vector{Float64}
    (p' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
end

"""
Apply affine transform to CartesianIndices
"""
function transform_point(tfm::AffineTransform,p::CartesianIndex{3})::Vector{Float64}
    ( [p[1]-1.0,p[2]-1.0,p[3]-1.0]' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
end


"""
Apply grid transform
"""
function transform_point(tfm::GridTransform, p::Vector{Float64})::Vector{Float64}
    # convert to voxel coordinates
    v=transform_point(tfm.world_to_voxel,p)
    # interpolate to the point
    # TODO
    g=[0.0,0.0,0.0]
    # final transform
    p+g
end

