# for Affine transforms
using LinearAlgebra

# for grid transforms
using Interpolations

# for quick operations
using StaticArrays

"""
Affine transform
"""
struct AffineTransform{T} 
    rot::SMatrix{3,3,T,9}
    shift::SVector{3,T}
    
    function AffineTransform{T}(rot::SMatrix{3,3,T,9},
                                shift::SVector{3,T}) where {T}
        new(rot, shift)
    end

    function AffineTransform{T}() where {T}
        new( SMatrix{3,3,T,9}( [1.0 0.0 0.0 ;0.0 1.0 0.0 ;0.0 0.0 1.0 ]), 
             SVector{3,T}( [0.0,0.0,0.0] ) )
    end

    function AffineTransform{T}(mat::Matrix{T}) where {T} 
        new( SMatrix{3,3,T,9}(mat[1:3,1:3]), SVector{3,T}(mat[1:3,4]))
    end

    function AffineTransform{T}(mat::SMatrix{4,4,T,16}) where {T} 
        new( SMatrix{3,3,T,9}(mat[1:3,1:3]), SVector{3,T}(mat[1:3,4]))
    end

    function AffineTransform{T}(rot::Matrix{T}, shift::Vector{T}) where {T} 
        new( SMatrix{3,3,T,9}(rot), SVector{3,T}(shift))
    end
end

"""
Dense vector field transform (grid transform)
"""
struct GridTransform{T,F}
    voxel_to_world::AffineTransform{T}
    world_to_voxel::AffineTransform{T}
    vector_field::Array{F, 4}
    itp_vector_field

    function GridTransform{T,F}(voxel_to_world::AffineTransform{T},
        vector_field::Array{F, 4}) where {T,F} 
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, 
                (NoInterp(),BSpline(Linear()),BSpline(Linear()),BSpline(Linear()))),Flat()))
    end
end

function GridTransform() where {T,F}
    GridTransform{T,F}(
        AffineTransform{T}(),
        zeros{F}(3,3,3,3)
    )
end

"""
Dense vector field transform (grid transform) used in inverse
"""
struct InverseGridTransform{T,F}
    voxel_to_world::AffineTransform{T}
    world_to_voxel::AffineTransform{T}
    vector_field::Array{F, 4}

    itp_vector_field

    function InverseGridTransform{T,F}(voxel_to_world::AffineTransform{T},
        vector_field::Array{F, 4}) where {T,F}
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, 
                (NoInterp(),BSpline(Linear()), BSpline(Linear()), BSpline(Linear()))),Flat()))
    end
end

function InverseGridTransform() where {T,F}
    InverseGridTransform{T,F}(
        AffineTransform{T}(),
        zeros{F}(3,3,3,3)
    )
end


"""
AnyTransform
"""
AnyTransform{T,F}=Union{AffineTransform{T}, GridTransform{T,F}, InverseGridTransform{T,F}}


"""
Invert AffineTransform transform
"""
function inv(t::AffineTransform{T})::AffineTransform{T} where {T,F}
    AffineTransform{T}(Base.inv( [t.rot t.shift;0 0 0 1] ))
end

"""
Invert GridTransform transform
"""
function inv(t::GridTransform{T,F})::InverseGridTransform{T,F} where {T,F}
    InverseGridTransform{T,F}(t.voxel_to_world, t.vector_field)
end

"""
Invert GridTransform transform
"""
function inv(t::InverseGridTransform{T,F})::GridTransform{T,F} where {T,F}
    GridTransform{T,F}(t.voxel_to_world, t.vector_field)
end


"""
Invert concatenated transform
"""
function inv(t::Vector{AnyTransform{T,F}})::Vector{AnyTransform{T,F}} where {T,F}
    [inv(i) for i in reverse(t)]
end

"""
Apply affine transform
"""
@inline function transform_point(tfm::AffineTransform{T}, 
        p::SVector{3,T}; 
        _whatever...)::SVector{3,T} where {T,F}
    
    (p' * tfm.rot)' + tfm.shift
end


@inline function interpolate_field(v2w::AffineTransform{T},
        itp_vector_field, p::SVector{3,T})::SVector{3,T} where {T,F}
    # convert to voxel coords, add 1 to get index
    v = transform_point{T}(v2w, p) + SVector{3,T}(1.0,1.0,1.0)
    SVector{3,T}(itp_vector_field(1,v...),
                 itp_vector_field(2,v...),
                 itp_vector_field(3,v...) )
end

"""
Apply forward grid transform
"""
@inline function transform_point(tfm::GridTransform{T,F}, p::SVector{3,T};
        max_iter::Int=10,ftol::Float64=1e-3)::SVector{3,T} where {T,F}
    return p + interpolate_field{T}(tfm.world_to_voxel, tfm.itp_vector_field)
end

"""
Apply inverse grid transform
reimplements algorithm from MNI_formats/grid_transforms.c:grid_inverse_transform_point
"""
@inline function transform_point(tfm::InverseGridTransform{T,F}, p::SVector{3,T};
        max_iter::Int=10,ftol::Float64=1.0/80)::SVector{3,T}  where {T,F}
    
    best::SVector{3,T} = estimate::SVector{3,T} = p-interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
    err::SVector{3,T} = p-(estimate+interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))

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
@inline function transform_point(tfm::Vector{AnyTransform{T,F}}, 
        p::SVector{3,T};
        max_iter::Int=10,ftol::Float64=1.0/80)::SVector{3,T} where {T,F}
    for t in tfm
        p = transform_point(t,p;max_iter,ftol)
    end
    return p
end


"""
Apply affine transform to CartesianIndices
"""
@inline function transform_point(tfm::AffineTransform{T}, p::CartesianIndex{3};
        _whatever...)::SVector{3,T} where {T,F}
    (SVector{3,T}(p[1]-1.0, p[2]-1.0, p[3]-1.0)' * tfm.rot)' + tfm.shift
end


# helper 
Base.show(io::IO, z::GridTransform) = print(io, "GridTransform:", size(z.vector_field))
# helper 
Base.show(io::IO, z::InverseGridTransform) = print(io, "InverseGridTransform:", size(z.vector_field))


# NOTES
# TODO: implement method from 
#  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2837727/ 
# "A pseudoinverse deformation vector field generator and its applications"
# or https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097910/
# "Iterative inversion of deformation vector fields with feedback control" 
