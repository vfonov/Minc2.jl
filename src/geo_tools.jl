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
    vector_field::Array{Float64, 4}
    itp_vector_field
    function GridTransform(voxel_to_world::AffineTransform,
        vector_field::Array{Float64, 4})
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, (NoInterp(),BSpline(Linear()),BSpline(Linear()),BSpline(Linear()))),Flat()))
    end
end

function GridTransform()
    GridTransform(
        AffineTransform(),
        zeros(3,3,3,3)
    )
end

"""
Dense vector field transform (grid transform) used in inverse
"""
mutable struct InverseGridTransform
    voxel_to_world::AffineTransform
    world_to_voxel::AffineTransform
    vector_field::Array{Float64, 4}
    itp_vector_field
    function InverseGridTransform(voxel_to_world::AffineTransform,
        vector_field::Array{Float64, 4})
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, (NoInterp(),BSpline(Linear()),BSpline(Linear()),BSpline(Linear()))),Flat()))
    end
end

function InverseGridTransform()
    InverseGridTransform(
        AffineTransform(),
        zeros(3,3,3,3)
    )
end


"""
AnyTransform
"""
AnyTransform=Union{AffineTransform,GridTransform,InverseGridTransform}


"""
Invert AffineTransform transform
"""
function inv(t::AffineTransform)::AffineTransform
    AffineTransform(Base.inv(t.mat))
end

"""
Invert GridTransform transform
"""
function inv(t::GridTransform)::InverseGridTransform
    InverseGridTransform(t.voxel_to_world, t.vector_field)
end

"""
Invert GridTransform transform
"""
function inv(t::InverseGridTransform)::GridTransform
    GridTransform(t.voxel_to_world, t.vector_field)
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


function interpolate_field(v2w::AffineTransform,
        itp_vector_field,p::Vector{Float64})::Vector{Float64}
    # convert to voxel coords, add 1 to get index
    v = transform_point(v2w, p) .+ 1.0
    [itp_vector_field(1,v...),
     itp_vector_field(2,v...),
     itp_vector_field(3,v...) ]
end

"""
Apply forward grid transform
"""
function transform_point(tfm::GridTransform, p::Vector{Float64})::Vector{Float64}
    return p+interpolate_field(tfm.world_to_voxel,tfm.itp_vector_field)
end

"""
Apply inverse grid transform
reimplements algorithm from MNI_formats/grid_transforms.c:grid_inverse_transform_point
"""
function transform_point(tfm::InverseGridTransform, p::Vector{Float64})::Vector{Float64}
    #(p' * tfm.mat[1:3, 1:3])' + tfm.mat[1:3,4]
    # convert to voxel coords, add 1 to get index
    ftol=1.0/80
    max_iter=10

    best = estimate = p-interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
    err = p-(estimate+interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))

    smallest_err=sum(abs.(err))
    i=1

    while i<max_iter && smallest_err>ftol 
        i+=1
        estimate = estimate + 0.95 * err
        err = p-(estimate+interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))
        err_mag=sum(abs.(err))

        if err_mag<smallest_err
            best=estimate
            err_mag<smallest_err
        end
    end

    return best
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
# helper 
Base.show(io::IO, z::InverseGridTransform) = print(io, "InverseGridTransform:", size(z.vector_field))
