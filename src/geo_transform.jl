# Geometric transformations functions

# for Affine transforms
using LinearAlgebra

# for grid transforms
using Interpolations

# for quick operations
using StaticArrays



"""
    AnyTransform()

Any 3D Geometrical transform, abstract type
"""
abstract type AnyTransform end

"""
    Base.length(tr::AnyTransform)

Length of concatenated transform
"""
Base.length(tr::AnyTransform)=1

"""
    IdentityTransform()

Identity transform, does nothing
"""
struct IdentityTransform <: AnyTransform
end


"""
    transform_point(tfm::IdentityTransform, 
        p::SVector{3,T};
        _whatever...)

* `tfm` - IdentityTransform
* `p` - point to transform

Apply identity transform to a point (returns same point)
"""
@inline function transform_point(tfm::IdentityTransform, 
    p::SVector{3,T};
    _whatever...)::SVector{3,T} where {T}
    p
end


"""
    AffineTransform{T}

Affine transform, described by rotation matrix and shift vector
"""
struct AffineTransform{T} <: AnyTransform
    rot::SMatrix{3,3,T,9}
    shift::SVector{3,T}
end
function AffineTransform(::Type{T}=Float64) where {T}
    return AffineTransform( SMatrix{3,3,T,9}( [1 0 0 ;0 1 0 ;0 0 1 ]), 
                             SVector{3,T}( [0,0,0] ) )
end

"""
    AffineTransform(mat)  Create affine transform from a matrix-like object

* `mat` - matrix  (4x4 or 3x4)
"""
function AffineTransform(mat) 
    ind = SA[1, 2, 3]
    return AffineTransform(mat[ind, ind], mat[ind, 4])
end

"""
    AffineTransform(rot, shift) Create affine transform from rotation and shift

* `rot` - rotation matrix-like (3x3)

* `shift` - shift vector-like (3)
"""
function AffineTransform(rot, shift)
     ind = SA[1, 2, 3]
     return AffineTransform(rot[ind, ind], shift[ind])
end


"""
    GridTransform{T,F,I} 

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

* `grid` - GridTransform

"""
voxel_to_world(grid::GridTransform) = grid.voxel_to_world


"""
    world_to_voxel(grid::GridTransform)

Extract world to voxel affine transform from a GridTransform

* `grid` - GridTransform
"""
world_to_voxel(grid::GridTransform) = grid.world_to_voxel

"""
    array(grid::GridTransform)
 
Extract underlying plain array

* `grid` - GridTransform
"""
array(grid::GridTransform) = grid.vector_field



"""
    GridTransform(
        voxel_to_world::AffineTransform{T},
        vector_field::Array{F, 4})

Construct GridTransform from voxel to world transform and a vector field

* `voxel_to_world` - voxel to world affine transform
* `vector_field` - vector field
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
    GridTransform(::Type{T}=Float64,::Type{F}=Float64)

Construct empty GridTransform, which should generate identity transform
"""
function GridTransform(::Type{T}=Float64,::Type{F}=Float64) where {T,F}
    GridTransform(
        AffineTransform(T),
        zeros(F,3,3,3,3)
    )
end

"""
    InverseGridTransform{T,F,VF}
    
Dense vector field transform (GridTransform) used in inverse
"""
struct InverseGridTransform{T,F,VF} <: AnyTransform
    voxel_to_world::AffineTransform{T}
    world_to_voxel::AffineTransform{T}
    vector_field::Array{F, 4}
    itp_vector_field::VF

end


"""
    voxel_to_world(grid::InverseGridTransform)

Extract voxel to world affine transform from a InverseGridTransform

* `grid` - InverseGridTransform
"""
voxel_to_world(grid::InverseGridTransform) = grid.voxel_to_world


"""
    world_to_voxel(grid::InverseGridTransform)

Extract world to voxel affine transform from a InverseGridTransform

* `grid` - InverseGridTransform
"""
world_to_voxel(grid::InverseGridTransform) = grid.world_to_voxel


"""
    array(grid::InverseGridTransform)

Extract underlying plain array

* `grid` - InverseGridTransform
"""
array(grid::InverseGridTransform) = grid.vector_field


"""
    InverseGridTransform(
        voxel_to_world::AffineTransform{T},
        vector_field::Array{F, 4})

Construct `InverseGridTransform` from voxel to world transform
    and a vector field

* `voxel_to_world` - voxel to world affine transform
* `vector_field` - vector field
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
    InverseGridTransform(::Type{T}=Float64,::Type{F}=Float64)

Construct `InverseGridTransform` empty transform
"""
function InverseGridTransform(::Type{T}=Float64,::Type{F}=Float64) where {T,F}
    InverseGridTransform(
        AffineTransform(T),
        zeros(F,3,3,3,3)
    )
end


"""
    GeoTransforms=Vector{AnyTransform}

Concatenated transforms
"""
GeoTransforms=Vector{AnyTransform}


"""
    inv(::IdentityTransform)::IdentityTransform

Invert IdentityTransform transform, does nothing
"""
function inv(::IdentityTransform)::IdentityTransform
    IdentityTransform()
end


"""
    inv(t::AffineTransform{T})::AffineTransform{T}

Invert AffineTransform transform

* `t` - AffineTransform
"""
function inv(t::AffineTransform{T})::AffineTransform{T} where {T}
    AffineTransform(Base.inv( SMatrix{4,4,T,16}([t.rot t.shift;0 0 0 1]) ))
end


"""
    inv(t::GridTransform{T,F,VF})::InverseGridTransform{T,F,VF}

Invert GridTransform transform

* `t` - GridTransform
"""
function inv(t::GridTransform{T,F,VF})::InverseGridTransform{T,F,VF} where {T,F,VF}
    InverseGridTransform( t.voxel_to_world, t.vector_field)
end


"""
    inv(t::InverseGridTransform{T,F,VF})::GridTransform{T,F,VF}

Invert InverseGridTransform transform

* `t` - InverseGridTransform
"""
function inv(t::InverseGridTransform{T,F,VF})::GridTransform{T,F,VF} where {T,F,VF}
    GridTransform(t.voxel_to_world, t.vector_field)
end


"""
    inv(t::Vector{T})::Vector{AnyTransform}

Invert concatenated transform

* `t` - concatenated transform
"""
function inv(t::Vector{T})::Vector{AnyTransform} where T<:AnyTransform
    AnyTransform[inv(i) for i in reverse(t)]
end


"""
    transform_point(
            tfm::AffineTransform{T}, 
            p::SVector{3,T};
            _whatever...)::SVector{3,T}

Apply affine transform to a point

* `tfm` - affine transform
* `p` - point to transform
"""
@inline function transform_point(
        tfm::AffineTransform{T}, 
        p::SVector{3,T};
        _whatever...)::SVector{3,T} where {T}
    
    return tfm.rot*p + tfm.shift
end


"""
    interpolate_field(
            v2w::AffineTransform{T},
            itp_vector_field::I, 
            p::SVector{3,T} )::SVector{3,T}

Internal support function

* `v2w` - voxel to world affine transform
* `itp_vector_field` - vector field
* `p` - point to transform
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
    transform_point(
            tfm::GridTransform{T,F}, p::SVector{3,T};
            _whatever...)::SVector{3,T} where {T,F}

Apply forward grid transform to a point

* `tfm` - grid transform
* `p` - point to transform
"""
@inline function transform_point(
        tfm::GridTransform{T,F}, p::SVector{3,T};
        _whatever...)::SVector{3,T} where {T,F}
    return p + interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
end


"""
    transform_point(
        tfm::AffineTransform{T}, 
        p::CartesianIndex{3};
        _whatever...)::SVector{3,T}

Apply inverse grid transform
reimplements algorithm from MNI_formats/grid_transforms.c:grid_inverse_transform_point

* `tfm` - affine transform
* `p` - point to transform
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
    transform_point(
        tfm::Vector{XFM},
        p::SVector{3,T};
        max_iter::Int=10,
        ftol::Float64=1.0/80)::SVector{3,T}

Apply concatenated transform to a point

* `tfm` - concatenated transform
* `p` - point to transform
* `max_iter` - maximum number of iterations for inverse transform
* `ftol` - tolerance for inverse transform
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
    transform_point(
        tfm::AffineTransform{T}, 
        p::CartesianIndex{3};
        _whatever...)::SVector{3,T}


Apply affine transform to CartesianIndices

* `tfm` - affine transform
* `p` - point to transform
"""
@inline function transform_point(
        tfm::AffineTransform{T}, 
        p::CartesianIndex{3};
        _whatever...)::SVector{3,T} where {T}
    tfm.rot * SVector{3,T}(p[1]-1.0, p[2]-1.0, p[3]-1.0) + tfm.shift
end


"""
   decompose(rot, shift)

Decompose affine transform specified as rotation matrix ans shift vector into three components
start, step, direction cosines

* `rot` - rotation matrix (3x3)
* `shift` - shift vector (3)
"""
function decompose(rot, shift)
    f = svd(rot)

    # remove scaling
    dir_cos = f.U * f.Vt

    step  = diag(Base.inv(dir_cos) * rot)
    start = Base.inv(dir_cos) * shift
    
    return start, step, dir_cos
end


"""
    decompose(tfm::AffineTransform{T})

Decompose affine transform into three components
start, step, direction cosines

* `tfm` - affine transform
"""
function decompose(tfm::AffineTransform{T}) where {T}
    return decompose(tfm.rot, tfm.shift)
end


"""
    decompose(tfm::Matrix{T})

Decompose affine transform into three components
start, step, direction cosines
"""
function decompose(tfm::Matrix{T}) where {T}
    return decompose( tfm[1:3, 1:3], tfm[1:3, 4])
end


# helper
"""
    Base.show(io::IO, z::GridTransform{T,F,I})

Print summary information about grid transform
"""
Base.show(io::IO, z::GridTransform{T,F,I}) where {T,F,I} = print(io, "GridTransform{$(T),$(F),...}:", size(z.vector_field))

# helper 
"""
    Base.show(io::IO, z::InverseGridTransform{T,F,I})

Print summary information about inverse grid transform
"""
Base.show(io::IO, z::InverseGridTransform{T,F,I})  where {T,F,I} = print(io, "InverseGridTransform{$(T),$(F),...}:", size(z.vector_field))


# NOTES
# TODO: implement method from 
#  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2837727/ 
# "A pseudoinverse deformation vector field generator and its applications"
# or https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097910/
# "Iterative inversion of deformation vector fields with feedback control" 
