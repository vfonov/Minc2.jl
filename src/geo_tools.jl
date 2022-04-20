"""
Affine transform
"""
mutable struct AffineTransform
    mat::Matrix{Float64}
    function AffineTransform()
        new([1.0 0.0 0.0 0.0;0.0 1.0 0.0 0.0;0.0 0.0 1.0 0.0;0.0 0.0 0.0 1.0])
    end
end


"""
Dense vector field transform (grid transform)
"""
mutable struct GridTransform
    world_to_voxel::AffineTransform
    voxel_to_world::AffineTransform
    inverse::Bool
    vector_field::Array{Float64, 4}
    function GridTransform()
        new(AffineTransform(),AffineTransform(),false,zeros(1,1,1,3))
    end
end

"""
AnyTransform
"""
AnyTransform=Union{AffineTransform,GridTransform}


"""
Apply affine transform
"""
function transform_point(tfm::AffineTransform,p::Vector{Float64})::Vector{Float64}
    (p' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
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

