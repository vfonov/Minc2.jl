# minimally working example for discourse.julialang.org/

using StaticArrays

struct AffineTransform{T} 
    rot::SMatrix{3,3,T,9}
    shift::SVector{3,T}
end

function AffineTransform(mat) 
    ind = SA[1, 2, 3]
    return AffineTransform(mat[ind, ind], mat[ind, 4])
end

function AffineTransform(rot, shift)
    ind = SA[1, 2, 3]
    return AffineTransform(rot[ind, ind], shift[ind])
end


# apply transform to a physical coordinate
@inline function transform_point(tfm::AffineTransform{T}, 
    p::SVector{3,T})::SVector{3,T} where {T}

    return (p' * tfm.rot)' + tfm.shift
end

# apply transform to an index, need to subtract 1 
@inline function transform_point(tfm::AffineTransform{T}, 
    p::CartesianIndex{3})::SVector{3,T} where {T}
    return ( SVector{3,T}(p[1]-1.0, p[2]-1.0, p[3]-1.0)' * tfm.rot)' + tfm.shift
end

#Invert AffineTransform transform
function inv(t::AffineTransform{T})::AffineTransform{T} where {T}
    AffineTransform(Base.inv( SMatrix{4,4,T,16}([t.rot t.shift;0 0 0 1]) ))
end


# actual example starts here
using Interpolations

size=256

# input volume
in_vol=rand(size,size,size)

# output volume
out_vol=similar(in_vol)

# transformation from o-based ontinious indices to physical coordinates:
i2w=AffineTransform(
    [1.0 0 0;0 1 0;0 0 1],[-size/2,-size/2,-size/2]
)

# transformation from physical coordinates back to continious indeces
w2i=inv(i2w)

# spatial transformation: rotate around axis by 30 degrees, no shift
tfm=AffineTransform(
    [
        1 0 0 0;
        0 0.866025388240814 -0.5 0;
        0 0.5 0.866025388240814 0;
        0 0 0 1
    ]
)

# setup intepolation object
itp = extrapolate(interpolate(in_vol, BSpline(Linear())), Flat())

# actual code that does computations:
function resample_volume(out_vol,i2w,w2i,itfm,itp)
    for c in CartesianIndices(out_vol)
        orig = transform_point(i2w, c )
        dst  = transform_point(itfm, orig)
        dst_v= transform_point(w2i, dst ) .+ 1.0
        
        @inbounds out_vol[c] = itp( dst_v... )
    end
    nothing
end

@timev resample_volume(out_vol,i2w,w2i,tfm,itp)
# TODO : run benchmarktools?
