# Geometric transformations functions

# for Affine transforms
using LinearAlgebra

# for grid transforms
using Interpolations

# for quick operations
using StaticArrays



"""
Any 3D Geometrical transform
"""
abstract type AnyTransform end


Base.length(tr::AnyTransform)=1

"""
Identity transform
"""
struct IdentityTransform <: AnyTransform
end

@inline function transform_point(tfm::IdentityTransform, 
    p::SVector{3,T};
    _whatever...)::SVector{3,T} where {T}
    p
end


"""
Affine transform, described by rotation matrix and shift vector
"""
struct AffineTransform{T} <: AnyTransform
    rot::SMatrix{3,3,T,9}
    shift::SVector{3,T}
end

"""
Create identity transform
"""
function AffineTransform(::Type{T}=Float64) where {T}
    return AffineTransform( SMatrix{3,3,T,9}( [1 0 0 ;0 1 0 ;0 0 1 ]), 
                             SVector{3,T}( [0,0,0] ) )
 end

"""
Create affine transform from a matrix-like object  (4x4 or 3x4)
"""
function AffineTransform(mat) 
    ind = SA[1, 2, 3]
    return AffineTransform(mat[ind, ind], mat[ind, 4])
end

"""
Create affine transform from rotation matrix-like object (3x3) and shift vector-like object (3)
"""
function AffineTransform(rot, shift)
     ind = SA[1, 2, 3]
     return AffineTransform(rot[ind, ind], shift[ind])
end


"""
Dense vector field transform (grid transform)
"""
struct GridTransform{T,F,I} <: AnyTransform
    voxel_to_world::AffineTransform{T}
    world_to_voxel::AffineTransform{T}
    vector_field::Array{F, 4}
    itp_vector_field::I
end


"""
    voxel_to_world(grid::GridTransform)

Extract voxel to world affine transform from a GridTransform
"""
voxel_to_world(grid::GridTransform) = grid.voxel_to_world


"""
    world_to_voxel(grid::GridTransform)

Extract world to voxel affine transform from a GridTransform
"""
world_to_voxel(grid::GridTransform) = grid.world_to_voxel

"""
    array(grid::GridTransform)
 
Extract underlying plain array
"""
array(grid::GridTransform) = grid.vector_field



"""
Construct GridTransform from voxel to world transform and a vector field
"""
function GridTransform(
    voxel_to_world::AffineTransform{T},
    vector_field::Array{F, 4}) where {T,F}

    GridTransform(voxel_to_world, inv(voxel_to_world), vector_field,
        extrapolate(interpolate(vector_field, 
                (NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()))),
            Flat()))
end


"""
Construct empty GridTransform, which should generate identity transform
"""
function GridTransform(::Type{T}=Float64,::Type{F}=Float64) where {T,F}
    GridTransform(
        AffineTransform(T),
        zeros(F,3,3,3,3)
    )
end

"""
Dense vector field transform (GridTransform) used in inverse
"""
struct InverseGridTransform{T,F,VF} <: AnyTransform
    voxel_to_world::AffineTransform{T}
    world_to_voxel::AffineTransform{T}
    vector_field::Array{F, 4}
    itp_vector_field::VF

end


"""
Extract voxel to world affine transform from a InverseGridTransform
"""
voxel_to_world(grid::InverseGridTransform) = grid.voxel_to_world


"""
Extract world to voxel affine transform from a InverseGridTransform
"""
world_to_voxel(grid::InverseGridTransform) = grid.world_to_voxel


"""
Extract underlying plain array
"""
array(grid::InverseGridTransform) = grid.vector_field


"""
Construct `InverseGridTransform` from voxel to world transform
    and a vector field
"""
function InverseGridTransform(
    voxel_to_world::AffineTransform{T},
    vector_field::Array{F, 4}) where {T,F}
    InverseGridTransform(voxel_to_world,inv(voxel_to_world),vector_field,
        extrapolate(interpolate(vector_field, 
                (NoInterp(), BSpline(Linear()), BSpline(Linear()), BSpline(Linear()))),
            Flat()))
end

"""
Construct `InverseGridTransform` empty transform
"""
function InverseGridTransform(::Type{T}=Float64,::Type{F}=Float64) where {T,F}
    InverseGridTransform(
        AffineTransform(T),
        zeros(F,3,3,3,3)
    )
end


"""
Concatenated transforms
"""
GeoTransforms=Vector{AnyTransform}


"""
Invert IdentityTransform transform
"""
function inv(::IdentityTransform)::IdentityTransform
    IdentityTransform()
end


"""
Invert AffineTransform transform
"""
function inv(t::AffineTransform{T})::AffineTransform{T} where {T}
    AffineTransform(Base.inv( SMatrix{4,4,T,16}([t.rot t.shift;0 0 0 1]) ))
end


"""
Invert GridTransform transform
"""
function inv(t::GridTransform{T,F,VF})::InverseGridTransform{T,F,VF} where {T,F,VF}
    InverseGridTransform( t.voxel_to_world, t.vector_field)
end


"""
Invert InverseGridTransform transform
"""
function inv(t::InverseGridTransform{T,F,VF})::GridTransform{T,F,VF} where {T,F,VF}
    GridTransform(t.voxel_to_world, t.vector_field)
end


"""
Invert concatenated transform
"""
function inv(t::Vector{T})::Vector{AnyTransform} where T<:AnyTransform
    AnyTransform[inv(i) for i in reverse(t)]
end


"""
Apply affine transform to a point
"""
@inline function transform_point(
        tfm::AffineTransform{T}, 
        p::SVector{3,T};
        _whatever...)::SVector{3,T} where {T}
    
    return (p' * tfm.rot)' + tfm.shift
end


"""
Internal support function
"""
@inline function interpolate_field(
        v2w::AffineTransform{T},
        itp_vector_field::I, 
        p::SVector{3,T} )::SVector{3,T} where {T, I<:Interpolations.Extrapolation}
    # convert to voxel coords, add 1 to get index
    v::SVector{3,T} = transform_point(v2w, p) .+ 1.0
    return SVector{3,T}(itp_vector_field(1,v...),
                        itp_vector_field(2,v...),
                        itp_vector_field(3,v...) )
end


"""
Apply forward grid transform to a point
"""
@inline function transform_point(
        tfm::GridTransform{T,F}, p::SVector{3,T};
        _whatever...)::SVector{3,T} where {T,F}
    return p + interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
end


"""
Apply inverse grid transform
reimplements algorithm from MNI_formats/grid_transforms.c:grid_inverse_transform_point
"""
@inline function transform_point(
        tfm::InverseGridTransform{T,F}, p::SVector{3,T};
        max_iter::Int=10, ftol::Float64=1.0/80)::SVector{3,T}  where {T,F}
    
    best::SVector{3,T} = estimate::SVector{3,T} = p - interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
    err::SVector{3,T} = p - (estimate + interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))

    smallest_err=sum(abs.(err))
    i=1

    while i<max_iter && smallest_err>ftol 
        i+=1
        estimate = estimate + 0.95 * err
        err = p - ( estimate + interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))
        err_mag=sum(abs.(err))

        if err_mag<smallest_err
            best = estimate
            err_mag<smallest_err
        end
    end

    return best
end


"""
Apply concatenated transform to a point
"""
@inline function transform_point(
        tfm::Vector{XFM},
        p::SVector{3,T};
        max_iter::Int=10,
        ftol::Float64=1.0/80)::SVector{3,T} where {XFM<:AnyTransform,T}
    for t in tfm
        p = transform_point(t,p;max_iter,ftol)
    end
    return p
end


# @inline function transform_point(
#     tfm::Tuple{},
#     p::SVector{3,T};
#     max_iter::Int=10,
#     ftol::Float64=1.0/80)::SVector{3,T} where {T}
#     return p
# end


# @inline function transform_point(
#     tfm::X,
#     p::SVector{3,T};
#     max_iter::Int=10,
#     ftol::Float64=1.0/80)::SVector{3,T} where {X<:Tuple,T}

#     transform_point(tfm[2:end],transform_point(tfm[1],p;max_iter,ftol);max_iter,ftol)
# end


"""
Apply affine transform to CartesianIndices
"""
@inline function transform_point(
        tfm::AffineTransform{T}, 
        p::CartesianIndex{3};
        _whatever...)::SVector{3,T} where {T}
    ( SVector{3,T}(p[1]-1.0, p[2]-1.0, p[3]-1.0)' * tfm.rot)' + tfm.shift
end


"""
Decompose affine transform into three components
start, step, direction cosines
"""
function decompose(rot, shift)
    f = svd(rot)

    # remove scaling
    dir_cos = f.U * f.Vt

    step  = diag(rot           * Base.inv(dir_cos))
    start = permutedims(permutedims(shift) * Base.inv(dir_cos))
    
    return start, step, dir_cos
end


"""
Decompose affine transform into three components
start, step, direction cosines
"""
function decompose(tfm::AffineTransform{T}) where {T}
    return decompose(tfm.rot, tfm.shift)
end


"""
Decompose affine transform into three components
start, step, direction cosines
"""
function decompose(tfm::Matrix{T}) where {T}
    return decompose(tfm[1:3,1:3], tfm[1:3,4])
end


# helper
"""
Print summary information about grid transform
"""
Base.show(io::IO, z::GridTransform{T,F,I}) where {T,F,I} = print(io, "GridTransform{$(T),$(F),...}:", size(z.vector_field))

# helper 
"""
Print summary information about grid transform
"""
Base.show(io::IO, z::InverseGridTransform{T,F,I})  where {T,F,I} = print(io, "InverseGridTransform{$(T),$(F),...}:", size(z.vector_field))


# NOTES
# TODO: implement method from 
#  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2837727/ 
# "A pseudoinverse deformation vector field generator and its applications"
# or https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097910/
# "Iterative inversion of deformation vector fields with feedback control" 
