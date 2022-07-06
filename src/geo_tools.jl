# for Affine transforms
using LinearAlgebra

# for grid transforms
using Interpolations

# for quick operations
using StaticArrays

"""
Affine transform
"""
struct AffineTransform
    rot::SMatrix{3,3,Float64}
    shift::SVector{3,Float64}
end

function AffineTransform()
    AffineTransform([1.0 0.0 0.0 ;0.0 1.0 0.0 ;0.0 0.0 1.0 ],[0.0,0.0,0.0])
end

function AffineTransform(mat::Matrix{Float64})
    AffineTransform(mat[1:3,1:3],mat[1:3,4])
end

function AffineTransform(mat::SMatrix{4,4,Float64})
    AffineTransform(mat[1:3,1:3],mat[1:3,4])
end


"""
Dense vector field transform (grid transform)
"""
struct GridTransform
    voxel_to_world::AffineTransform
    world_to_voxel::AffineTransform
    vector_field::Array{Float64, 4}
    itp_vector_field::Interpolations.Extrapolation{Float64, 4, 
        Interpolations.BSplineInterpolation{Float64, 4, Array{Float64, 4}, 
        Tuple{NoInterp, BSpline{Linear{Throw{OnGrid}}}, 
                        BSpline{Linear{Throw{OnGrid}}}, 
                        BSpline{Linear{Throw{OnGrid}}}}, 
        NTuple{4, Base.OneTo{Int64}}}, Tuple{NoInterp, BSpline{Linear{Throw{OnGrid}}}, 
        BSpline{Linear{Throw{OnGrid}}}, BSpline{Linear{Throw{OnGrid}}}}, Flat{Nothing}}

    function GridTransform(voxel_to_world::AffineTransform,
        vector_field::Array{Float64, 4})
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, 
                (NoInterp(),BSpline(Linear()),BSpline(Linear()),BSpline(Linear()))),Flat()))
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
struct InverseGridTransform
    voxel_to_world::AffineTransform
    world_to_voxel::AffineTransform
    vector_field::Array{Float64, 4}

    itp_vector_field::Interpolations.Extrapolation{Float64, 4, 
        Interpolations.BSplineInterpolation{Float64, 4, Array{Float64, 4}, 
        Tuple{NoInterp, BSpline{Linear{Throw{OnGrid}}}, 
                        BSpline{Linear{Throw{OnGrid}}}, 
                        BSpline{Linear{Throw{OnGrid}}}}, 
        NTuple{4, Base.OneTo{Int64}}}, Tuple{NoInterp, BSpline{Linear{Throw{OnGrid}}}, 
        BSpline{Linear{Throw{OnGrid}}}, BSpline{Linear{Throw{OnGrid}}}}, Flat{Nothing}}

    function InverseGridTransform(voxel_to_world::AffineTransform,
        vector_field::Array{Float64, 4})
        new(voxel_to_world,inv(voxel_to_world),vector_field,
            extrapolate(interpolate(vector_field, 
                (NoInterp(),BSpline(Linear()), BSpline(Linear()), BSpline(Linear()))),Flat()))
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
AnyTransform=Union{AffineTransform, GridTransform, InverseGridTransform}


"""
Invert AffineTransform transform
"""
function inv(t::AffineTransform)::AffineTransform
    AffineTransform(Base.inv( [t.rot t.shift;0 0 0 1] ))
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
@inline function transform_point(tfm::AffineTransform, 
        p::SVector{3,Float64}; 
        _whatever...)::SVector{3,Float64}
    
    (p' * tfm.rot)' + tfm.shift
end


@inline function interpolate_field(v2w::AffineTransform,
        itp_vector_field,p::SVector{3,Float64})::SVector{3,Float64}
    # convert to voxel coords, add 1 to get index
    v = transform_point(v2w, p) + SA_F64[1.0,1.0,1.0]
    SA_F64[itp_vector_field(1,v...),
           itp_vector_field(2,v...),
           itp_vector_field(3,v...) ]
end

"""
Apply forward grid transform
"""
@inline function transform_point(tfm::GridTransform, p::SVector{3,Float64};
        max_iter::Int=10,ftol::Float64=1e-3)::SVector{3,Float64}
    return p + interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field)
end

"""
Apply inverse grid transform
reimplements algorithm from MNI_formats/grid_transforms.c:grid_inverse_transform_point
"""
@inline function transform_point(tfm::InverseGridTransform, p::SVector{3,Float64};
        max_iter::Int=10,ftol::Float64=1.0/80)::SVector{3,Float64}
    
    best::SVector{3,Float64} = estimate::SVector{3,Float64} = p-interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, p)
    err::SVector{3,Float64} = p-(estimate+interpolate_field(tfm.world_to_voxel, tfm.itp_vector_field, estimate))

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
@inline function transform_point(tfm::Vector{AnyTransform}, p::SVector{3,Float64};
        max_iter::Int=10,ftol::Float64=1.0/80)::SVector{3,Float64}
    for t in tfm
        p::SVector{3,Float64} = transform_point(t,p;max_iter,ftol)
    end
    return p
end


"""
Apply affine transform to CartesianIndices
"""
@inline function transform_point(tfm::AffineTransform, p::CartesianIndex{3};
        _whatever...)::SVector{3,Float64}
    ( SA_F64[p[1]-1.0, p[2]-1.0, p[3]-1.0]' * tfm.rot)' + tfm.shift
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
