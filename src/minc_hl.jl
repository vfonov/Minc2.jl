# for resampling
using Interpolations
# for minc history
using Dates

"""
    Volume3D{T,N}

An abstract 3D volume, could be vector field or volume with time dimension
"""
struct Volume3D{T,N}
    "underlying array"
    vol::Array{T,N} 
    "Voxel to world affine transform"
    v2w::AffineTransform{Float64}
    "minc history"
    history::Union{String,Nothing} 
end

"""
    voxel_to_world(vol::Volume3D)

Extract voxel to world affine transform from a Volume3D
"""
voxel_to_world(vol::Volume3D) = vol.v2w

"""
    world_to_voxel(vol::Volume3D)

Extract world to voxel affine transform from a Volume3D
"""
function world_to_voxel(vol::Volume3D)
    v2w=voxel_to_world(vol)
    return inv(v2w)
end


"""
    array(vol::Volume3D)

Extract underlying plain array
"""
array(vol::Volume3D) = vol.vol


"""
    history(vol::Volume3D)

Extract history metadata
"""
history(vol::Volume3D) = vol.history


"""
    Volume3D(vol::Array{T,N}, v2w::AffineTransform{Float64}; 
        history::Union{AbstractString,Nothing}=nothing)::Volume3D{T,N}

Create Volume3D from an array and affine transform

* `vol` - underlying array
* `v2w` - voxel to world affine transform
* `history` - minc history
"""
function Volume3D(vol::Array{T,N}, v2w::AffineTransform{Float64}; 
        history::Union{AbstractString,Nothing}=nothing)::Volume3D{T,N} where {T,N}
    return Volume3D(vol, v2w, history)
end


"""
    Volume3D(vol::Array{T,N}, like::Volume3D; 
        history::Union{AbstractString,Nothing}=nothing)::Volume3D{T,N}

Create Volume3D from an array and another Volume3D that is used for sampling information

* `vol` - underlying array
* `like` - Volume3D that is used for sampling information
* `history` - minc history
"""
function Volume3D(vol::Array{T,N}, like::Volume3D; 
        history::Union{AbstractString,Nothing}=nothing)::Volume3D{T,N}  where {T,N}
    return Volume3D(vol, like.v2w, isnothing(history) ? like.history : history)
end


"""
    GridTransform(vol::V)

Create GridTransform from Volume3D, assume it is a vector field

* `vol` - Volume3D with 4D array containing vector field
"""
function GridTransform(vol::V) where {V<:Volume3D}
    return GridTransform(vol.v2w, vol.vol)
end


"""
    InverseGridTransform(vol::V)

Create InverseGridTransform from Volume3D, assume it is a vector field

* `vol` - Volume3D with 4D array containing vector field
"""
function InverseGridTransform(vol::V) where {V<:Volume3D}
    return InverseGridTransform(vol.v2w, vol.vol)
end


"""
    Base.show(io::IO, z::Volume3D{T,N})

Print Volume3D info
"""
Base.show(io::IO, z::Volume3D{T,N}) where {T,N} = print(io, "Volume3D{$(T),$(N)}:", size(z.vol)," ", z.v2w)


"""
    read_volume(fn::String; store::Type{T}=Float64)::Volume3D{T}

Read Volume3D from minc file

* `fn` - filename
* `store` - underlying array type
"""
function read_volume(fn::String; store::Type{T}=Float64)::Volume3D{T} where {T}
    in_vol,in_hdr,in_store_hdr,in_history = read_minc_volume_std_history(fn, store)
    v2w = voxel_to_world(in_hdr)

    return Volume3D(in_vol, v2w, in_history)
end


"""
    empty_volume_like(
            fn::String; 
            store::Type{T}=Float64, history=nothing)::Volume3D{T}

Create an empty Volume3D

* `fn` - filename of a minc file that is used for sampling information
* `store` - underlying array type
"""
function empty_volume_like(
        fn::String; 
        store::Type{T}=Float64, history=nothing)::Volume3D{T} where {T}
    out_vol,out_hdr, out_store_hdr, ref_history = empty_like_minc_volume_std_history(fn,store)
    v2w=voxel_to_world(out_hdr)

    return Volume3D(out_vol, v2w, isnothing(history) ? ref_history : history)
end


"""
    empty_volume_like(
        vol::Volume3D{T1,N}; 
        store::Type{T}=Float64, 
        history=nothing)

Create an empty Volume3D

* `vol` - Volume3D that is used for sampling information
* `store` - underlying array type
* `history` - minc history
"""
function empty_volume_like(
        vol::Volume3D{T1,N}; 
        store::Type{T}=Float64, 
        history=nothing) where {T1,T,N}
    out_vol = similar(vol.vol, store)
    return Volume3D(out_vol, vol.v2w, isnothing(history) ? vol.history : history )
end


"""
    full_volume_like(
            fn::String,
            fill::T=zero(T); 
            store::Type{T}=Float64, history=nothing)::Volume3D{T}

Create an empty Volume3D

* `fn` - filename of a minc file that is used for sampling information
* `store` - underlying array type
"""
function full_volume_like(
        fn::String, fill::T=zero(T);
        store::Type{T}=Float64, history=nothing)::Volume3D{T} where {T}
    out_vol,out_hdr, out_store_hdr, ref_history = empty_like_minc_volume_std_history(fn,store)
    v2w=voxel_to_world(out_hdr)
    v=Volume3D(out_vol, v2w, isnothing(history) ? ref_history : history)
    v.vol .= fill
    return v
end


"""
    empty_volume_like(
        vol::Volume3D{T1,N}; 
        store::Type{T}=Float64, 
        history=nothing)

Create an empty Volume3D

* `vol` - Volume3D that is used for sampling information
* `store` - underlying array type
* `history` - minc history
"""
function full_volume_like(
        vol::Volume3D{T1,N},fill::T=zero(T); 
        store::Type{T}=Float64, 
        history=nothing) where {T1,T,N}
    out_vol = similar(vol.vol, store)
    v= Volume3D(out_vol, vol.v2w, isnothing(history) ? vol.history : history )
    v.vol .= fill
    return v
end



"""
    save_volume(fn::AbstractString, 
        vol::Volume3D{T,N}; 
        store::Type{S}=Float32,
        history=nothing)

Save Volume3D to minc file

* `fn` - filename
* `vol` - Volume3D to save
* `store` - underlying MINC data type, to be used for storage
"""
function save_volume(fn::AbstractString, 
        vol::Volume3D{T,N}; 
        store::Type{S}=Float32,
        history=nothing) where {S,T,N}
    
    if isnothing(history)
        _history=vol.history
    else
        _history=(isnothing(vol.history) ? "" : vol.history * "\n" )*history
    end

    write_minc_volume_std(fn, store, 
      create_header_from_v2w(size(vol.vol), vol.v2w, 
      vector_dim=(length(size(vol.vol))==4)), vol.vol; history=_history)
end


"""
    resample_grid_volume!(
        in_vol::AbstractArray{T,4},
        out_vol::AbstractArray{T,4},
        v2w::AffineTransform{C}, 
        w2v::AffineTransform{C}, 
        itfm::Union{Vector{XFM}, XFM};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        fill=0.0,
        ftol=1.0/80,
        max_iter=10)::AbstractArray{T,4}

Resample 4D array using transformation,
assume 1st dimension is non spatial (vector dimension)

* `in_vol` - input 4D array
* `out_vol` - output 4D array
* `v2w` - voxel to world affine transform in the output array
* `w2v` - world to voxel affine transform in the input array
* `itfm` - inverse of the  transformation to apply (i.e from output to input)
* `interp` - interpolation method
* `fill` - fill value
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function resample_grid_volume!(
        in_vol::AbstractArray{T,4},
        out_vol::AbstractArray{T,4},
        v2w::AffineTransform{C}, 
        w2v::AffineTransform{C}, 
        itfm::Union{Vector{XFM}, XFM};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        fill=0.0,
        ftol=1.0/80,
        max_iter=10)::AbstractArray{T,4} where {T, C, I, XFM<:AnyTransform}

    # NEED to interpolate only over spatial dimensions
    in_vol_itp = extrapolate( interpolate( in_vol, (NoInterp(), interp, interp, interp)), fill)

    @simd for c in CartesianIndices(size(out_vol)[2:end])
        orig = transform_point(v2w, c )
        dst  = transform_point(itfm, orig; ftol, max_iter )
        dst_v= transform_point(w2v, dst ) .+ 1.0

        for i in eachindex(view(out_vol,:,1,1,1))
            @inbounds out_vol[i,c] = in_vol_itp(i, dst_v...)
        end
    end
    out_vol
end


"""
    resample_grid(
        in_grid::Volume3D{T,4}, 
        itfm::Union{Vector{XFM}, XFM}; 
        like::Union{Nothing,Volume3D{L,4}}=nothing)::Volume3D{T,4}

Resample Volume3D that contain 4D array,
using transformation, assume 1st dimension is non spatial

* `in_grid` - input Volume3D with 4D array describing vector field
* `itfm` - inverse of the  transformation to apply (i.e from output to input)
* `like` - Volume3D that is used for sampling information
"""
function resample_grid(
        in_grid::Volume3D{T,4}, 
        itfm::Union{Vector{XFM}, XFM}; 
        like::Union{Nothing,Volume3D{L,4}}=nothing)::Volume3D{T,4} where {T,L, XFM<:AnyTransform}

    if isnothing(like)
      out_vol = similar(in_grid.vol)
      v2w = in_grid.v2w
    else
      out_vol = similar(like.vol)
      v2w = like.v2w
    end
    resample_grid_volume!(in_grid.vol, out_vol, in_grid.v2w, inv(v2w), itfm; 
        interp=BSpline(Linear()))
    return Volume3D(out_vol, v2w)
end


"""
    tfm_to_grid!(
        tfm::Union{Vector{XFM}, XFM}, 
        grid::AbstractArray{T,4},
        v2w::AffineTransform{C};
        ftol=1.0/80,max_iter=10)::AbstractArray{T,4}

Convert arbitrary transformation 
into vector field

* `tfm` - transformation to use
* `grid` - output 4D array, will contain vector field
* `v2w` - voxel to world affine transform in the output array
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function tfm_to_grid!(
        tfm::Union{Vector{XFM}, XFM}, 
        grid::AbstractArray{T,4},
        v2w::AffineTransform{C};
        ftol=1.0/80,max_iter=10)::AbstractArray{T,4} where {T, C, XFM<:AnyTransform}
    
    @simd for c in CartesianIndices(size(grid)[2:end])
        orig = transform_point(v2w, c )
        dst  = transform_point(tfm, orig;max_iter,ftol)
        @inbounds grid[:,c] .= dst .- orig
    end
    return grid
end


"""
    tfm_to_grid(
        tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::Volume3D{T,4}

Convert arbitrary transformation 
into Volume3D with 4D array containing vector field

* `tfm` - transformation to use
* `v2w` - voxel to world affine transform in the output array
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function tfm_to_grid(
        tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::Volume3D{T,4} where {T, XFM<:AnyTransform, G<:Volume3D}
    # TODO: deal with 3D ref ?
    if ndims(ref.vol)==3
        # need to generate 4D output
        out_grid = Array{T,4}(undef,(3,size(ref.vol)...))
    else
        out_grid = similar(ref.vol, store)
    end
    v2w = ref.v2w

    tfm_to_grid!(tfm, out_grid, v2w; ftol, max_iter)
    return Volume3D( out_grid, v2w)
end


"""
    normalize_tfm(tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::GridTransform{Float64,T}

Convert arbitrary transformation  into a single GridTransform

* `tfm` - transformation to use
* `ref` - GridTransform that is used for sampling information
* `store` - underlying array type
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function normalize_tfm(tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::GridTransform{Float64,T} where {T, XFM<:AnyTransform, G<:GridTransform}

    out_grid = similar(ref.vector_field, store)
    v2w = voxel_to_world(ref)

    tfm_to_grid!(tfm,out_grid,v2w;ftol,max_iter)

    return GridTransform(v2w, out_grid)
end


"""
    normalize_tfm(tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::GridTransform{Float64,T}

Convert arbitrary transformation  into a single GridTransform

* `tfm` - transformation to use
* `ref` - Volume3D that is used for sampling information
* `store` - underlying array type
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function normalize_tfm(tfm::Union{Vector{XFM}, XFM},
        ref::G;
        store::Type{T}=Float64,ftol=1.0/80,max_iter=10)::GridTransform{Float64,T} where {T, XFM<:AnyTransform, G<:Volume3D}

    out_grid = similar(ref.vol, store)
    v2w = voxel_to_world(ref)

    tfm_to_grid!(tfm, out_grid, v2w;ftol,max_iter)

    return GridTransform(v2w, out_grid)
end


"""
    resample_volume!(in_vol::AbstractArray{T,3}, 
        out_vol::AbstractArray{T,3}, 
        v2w::AffineTransform{C}, 
        w2v::AffineTransform{C}, 
        itfm::Union{Vector{XFM},XFM};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        fill=0.0,
        ftol=1.0/80,
        max_iter=10)

Resample 3D array using transformation 

* `in_vol` - input 3D array
* `out_vol` - output 3D array
* `v2w` - voxel to world affine transform in the output array
* `w2v` - world to voxel affine transform in the input array
* `itfm` - inverse of the  transformation to apply (i.e from output to input)
* `interp` - interpolation method
* `fill` - fill value
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function resample_volume!(in_vol::AbstractArray{T,3}, 
        out_vol::AbstractArray{T,3}, 
        v2w::AffineTransform{C}, 
        w2v::AffineTransform{C}, 
        itfm::Union{Vector{XFM},XFM};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        fill=0.0,
        ftol=1.0/80,
        max_iter=10) where {C, T, I, XFM<:AnyTransform}

    in_vol_itp = extrapolate( interpolate( in_vol, interp), fill)
    @simd for c in CartesianIndices(out_vol)
        orig = transform_point(v2w, c )
        dst  = transform_point(itfm, orig; ftol, max_iter )
        dst_v= transform_point(w2v, dst ) .+ 1.0

        @inbounds out_vol[c] = in_vol_itp( dst_v... )
    end
    
    return out_vol
end


"""
    resample_volume!(
        in_vol::Volume3D{T,3}, 
        out_vol::Volume3D{O,3}; 
        tfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        itfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        interp::I=nothing, 
        fill=0.0, 
        order=nothing,
        ftol=1.0/80,
        max_iter=10)::Volume3D{O,3}

Resample Volume3D using transformation 
* `in_vol` - input Volume3D
* `out_vol` - output Volume3D
* `itfm` - inverse of the  transformation to apply (i.e from output to input)
* `tfm`  - transformation to apply (i.e from output to input) (instead of `itfm`)
* `interp` - interpolation method
* `fill` - fill value
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function resample_volume!(
        in_vol::Volume3D{T,3}, 
        out_vol::Volume3D{O,3}; 
        tfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        itfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        interp::I=nothing, 
        fill=0.0, 
        order=nothing,
        ftol=1.0/80,
        max_iter=10)::Volume3D{O,3} where {T,O,I,XFM<:AnyTransform}

    @assert ndims(out_vol.vol)==3
    @assert ndims(in_vol.vol)==3

    if !isnothing(tfm) && isnothing(itfm)
        itfm=inv(tfm)
    elseif isnothing(itfm) && isnothing(tfm)
        itfm=IdentityTransform() # identity transform
    end

    if isnothing(interp)
        if O <: Integer # output is integer, use nearest neighbor
            interp=BSpline(Constant())
            if !isnothing(order) && order!=0
                @error "Unsupported order when interpolating integers" order
            end
        else 
            if order==1
                interp=BSpline(Linear())
            elseif isnothing(order) || order==2  # default
                interp=BSpline(Quadratic(Line(OnCell())))
            elseif order==3
                interp=BSpline(Cubic(Line(OnCell())))
            else
                @error "Unsupported order" order
            end
        end
    end

    resample_volume!(in_vol.vol, out_vol.vol, out_vol.v2w, inv(in_vol.v2w), itfm; 
        interp, ftol, max_iter, fill)

    return out_vol
end


"""
    resample_volume(
        in_vol::Volume3D{T,3};
        like::Union{Volume3D{O,3},Nothing}=nothing,
        tfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        itfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        interp::I=nothing, 
        fill=0.0,
        order=1,
        ftol=1.0/80,
        max_iter=10)::Volume3D 

Resample Volume3D using transformation 

* `in_vol` - input Volume3D
* `like` - Volume3D that is used for sampling information
* `itfm` - inverse of the  transformation to apply (i.e from output to input)
* `tfm`  - transformation to apply (i.e from output to input) (instead of `itfm`)
* `interp` - interpolation method
* `fill` - fill value
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function resample_volume(
        in_vol::Volume3D{T,3};
        like::Union{Volume3D{O,3},Nothing}=nothing,
        tfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        itfm::Union{Vector{XFM},XFM,Nothing}=nothing, 
        interp::I=nothing, 
        fill=0.0,
        order=1,
        ftol=1.0/80,
        max_iter=10)::Volume3D where {T, O, I, XFM<:AnyTransform}

    out_vol = empty_volume_like( isnothing(like) ? in_vol : like ;
        store = (isnothing(like) ? eltype(in_vol.vol) : eltype(like.vol) ) )

    return resample_volume!(in_vol, out_vol; tfm, itfm, interp, fill, order, ftol, max_iter)
end


"""
    crop_volume(in_vol::Volume3D{T,N},crop;
        fill_val::T=zero(T))::Volume3D{T,N}

Crop (or pad) a volume

* `in_vol` - input Volume3D
* `crop` - crop specification, e.g. `[(1,2),(3,4),(5,6)]` , negative values mean padding
"""
function crop_volume(in_vol::Volume3D{T,N},crop;
        fill_val::T=zero(T))::Volume3D{T,N} where {T,N}

    sz = size(in_vol.vol)
    if N==4 # we have vector dimension (?)
        shift=1
    else
        shift=0
    end

    @assert (ndims(in_vol.vol)-shift) == length(crop) "Unequal dimrange size"

    new_sz    = [ sz[j+shift]-i[1]-i[2] for (j,i) in enumerate(crop)]

    in_range  = [(max(i[1]+1,1) : min(sz[j+shift]-i[2], sz[j+shift])) for (j,i) in enumerate(crop)]
    out_range = [(max(-i[1]+1,1): min(new_sz[j]+i[2],     new_sz[j])) for (j,i) in enumerate(crop)]

    if shift == 1
        new_sz    = [1;new_sz]
        in_range  = [1:sz[1]; in_range]
        out_range = [1:sz[1]; out_range]
    end

    old_v2w = voxel_to_world(in_vol)
    #start,step,dir_cos = decompose(in_vol.v2w)
    new_start = transform_point(old_v2w, SVector{3,Float64}( [crop[i][1] for i in 1:(N-shift)]))
    old_start = transform_point(old_v2w, SVector{3,Float64}( [0.0,0.0,0.0] ) )

    coord_shift = new_start - old_start
    out_array = fill(fill_val, new_sz...)

    @show new_sz out_range in_range

    out_array[out_range...] .= in_vol.vol[in_range...]

    return Volume3D(out_array, AffineTransform(old_v2w.rot, old_v2w.shift+coord_shift))
end


"""
    calculate_jacobian!(
        tfm::Union{Vector{XFM},XFM},
        out_vol::AbstractArray{T,3},
        out_v2w::AffineTransform{C};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10)

Calculate dense jacobian determinant field for an arbitrary transformation

* `tfm`  - transformation to use
* `out_vol` - output 3D array
* `out_v2w` - voxel to world affine transform in the output array
* `interp` - interpolation method
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function calculate_jacobian!(
        tfm::Union{Vector{XFM},XFM},
        out_vol::AbstractArray{T,3},
        out_v2w::AffineTransform{C};
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10) where {C,T,I,XFM<:AnyTransform}

    # calculate scaling matrix from the voxel to world matrix
    # to compensate for the step size (makes no difference on 1x1x1 voxels)
    _,step,_ = decompose(out_v2w)
    sc = diagm(inv.(step)) 

    # First step: generate vector field of transformations
    vector_field = Array{T}(undef, 3, size(out_vol)...)
    @inbounds let 
        @simd for c in CartesianIndices(out_vol)
            orig = transform_point(out_v2w, c )
            dst  = transform_point(tfm, orig; ftol, max_iter )

            vector_field[:,c] .= dst # .- orig # mincblob convention (?)
        end
    end
    # Second step: calculate jacobian determinant 
    vector_field_itp = extrapolate( interpolate( vector_field, 
        (NoInterp(),interp,interp,interp)), Flat())

    @simd for c in CartesianIndices(out_vol)
        grad = hcat([ Interpolations.gradient( vector_field_itp, i, Tuple(c)...) for i in 1:3 ]...)
        @inbounds out_vol[c] = det(grad'*sc)
    end

    return out_vol
end

"""
    calculate_jacobian!(
        tfm::Union{Vector{XFM},XFM}, 
        out_vol::Volume3D{T,3}; 
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10)::Volume3D{T,3}

Calculate dense jacobian determinant field for an arbitrary transformation

* `tfm`  - transformation to use
* `out_vol` - output Volume3D
* `interp` - interpolation method
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function calculate_jacobian!(
        tfm::Union{Vector{XFM},XFM}, 
        out_vol::Volume3D{T,3}; 
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10)::Volume3D{T,3} where {T,I,XFM<:AnyTransform}

    calculate_jacobian!(tfm, out_vol.vol,out_vol.v2w;interp,ftol,max_iter) 
    return out_vol
end


"""
    calculate_jacobian(
        tfm::Union{Vector{XFM},XFM};
        like::Volume3D{T,3}; 
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10)::Volume3D{T,3}

Calculate dense jacobian determinant field for an arbitrary transformation

* `tfm`  - transformation to use
* `out_vol` - output Volume3D
* `interp` - interpolation method
* `ftol` - tolerance, for inverse transformations
* `max_iter` - maximum number of iterations, for inverse transformations
"""
function calculate_jacobian(
        tfm::Union{Vector{XFM},XFM}, 
        like::Volume3D{T,3}; 
        interp::I=BSpline(Quadratic(Line(OnCell()))),
        ftol=1.0/80,
        max_iter=10)::Volume3D{T,3} where {T,I,XFM<:AnyTransform}

    out_vol=empty_volume_like(like)
    calculate_jacobian!(tfm, out_vol.vol,out_vol.v2w;interp,ftol,max_iter) 
    return out_vol
end



"""
    format_history(args)::String

Generate minc-style history from program args
"""
function format_history(args)::String
    return Dates.format(now(),"e d u HH:MM:SS YYYY")*">>>"*join(args," ")
end

