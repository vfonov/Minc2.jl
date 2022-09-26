# for resampling
using Interpolations


"""
An abstract 3D volume
"""
struct Volume3D
    vol # an abstract volume
    v2w # transformation from voxel to world coordinates
end

function read_volume(fn::String; store::Type{T}=Float64) where {T}
    in_vol,in_hdr,in_store_hdr = Minc2.read_minc_volume_std(fn, store)
    v2w=Minc2.voxel_to_world(in_hdr)

    return Volume3D(in_vol,v2w)
end


function empty_volume_like(fn::String; store::Type{T}=Float64) where {T}
    out_vol,out_hdr,out_store_hdr = Minc2.empty_like_minc_volume_std(fn,store)
    v2w=Minc2.voxel_to_world(out_hdr)

    return Volume3D(out_vol,v2w)
end


function empty_volume_like(vol::Volume3D; store::Type{T}=Float64) where {T}
    out_vol = Array{Float64}(undef, size(vol.vol)...)
    return Volume3D(out_vol,vol.v2w)
end


function save_volume(fn, vol::Volume3D; store::Type{T}=Float32) where {T}
    Minc2.write_minc_volume_std(fn, store, 
      Minc2.create_header_from_v2w(size(vol.vol), vol.v2w, vector_dim=(length(size(vol.vol))==4)), vol.vol)
end



function resample_volume(in_vol::Array{T,3}, 
    out_vol::Array{T,3}, 
    v2w::Minc2.AffineTransform{C}, 
    w2v::Minc2.AffineTransform{C}, 
    itfm::Vector{Minc2.AnyTransform};
    interp::I=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10) where {C,T,I}

    in_vol_itp = extrapolate( interpolate( in_vol, interp),fill)

    @simd for c in CartesianIndices(out_vol)
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(itfm, orig; ftol, max_iter )
        dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0

        @inbounds out_vol[c] = in_vol_itp( dst_v... )
        #out_vol[c] = sqrt(sum((orig - dst).^2))
    end
    return out_vol
end



# TODO merge with next
function resample_grid_volume(
    in_vol::Array{T,4}, 
    out_vol::Array{T,4}, 
    v2w::Minc2.AffineTransform{C}, 
    w2v::Minc2.AffineTransform{C}, 
    itfm::Vector{Minc2.AnyTransform};
    interp::I=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10) where {C,T,I}

    # NEED to interpolate only over spatial dimensions
    in_vol_itp = extrapolate( interpolate( in_vol, (NoInterp(),interp,interp,interp)),fill)

    @simd for c in CartesianIndices(view(out_vol,1,:,:,:))
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(itfm, orig; ftol, max_iter )
        dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0

        for i in eachindex(view(out_vol,:,1,1,1))
            @inbounds out_vol[i,c] = in_vol_itp(i, dst_v...)
        end
    end
    out_vol
end


# TODO: add function to apply local jacobian?
function resample_grid(in_grid, itfm; ref=nothing)
    if isnothing(ref)
      out_vol = similar(in_grid.vol)
      v2w = in_grid.v2w
    else
      out_vol = similar(ref.vol)
      v2w = ref.v2w
    end
  
    interp = BSpline(Linear())
    in_vol_itp = extrapolate( interpolate( in_grid.vol, (NoInterp(), interp, interp, interp)), 0.0)
    w2v=Minc2.inv(v2w)
  
    @simd for c in CartesianIndices(view(out_vol,1,:,:,:))
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(itfm, orig)
        dst_v= Minc2.transform_point(w2v, dst ) .+ 1.0
        
        for i in eachindex(view(out_vol,:,1,1,1))
            @inbounds out_vol[i,c] = in_vol_itp(i, dst_v...)
        end
    end
    return Volume3D(vol=out_vol, v2w=v2w)
end
  
# convert transforms into a single nonlinear grid transform
function normalize_tfm(tfm, ref)
    out_grid = similar(ref.vol)
    v2w = ref.v2w
  
    @simd for c in CartesianIndices(view(out_grid,1,:,:,:))
      orig = Minc2.transform_point(v2w, c )
      dst  = Minc2.transform_point(tfm, orig)
  
      out_grid[:,c] .= dst .- orig
    end
  
    return Minc2.GridTransform{Float64,Float64}(v2w, out_grid)
end


#TODO: make this more generic , and make it work with labels 
function resample_volume!(out_vol::Volume3D, in_vol::Volume3D; 
    tfm=nothing, itfm=nothing, interp=nothing, fill=nothing, 
    order=nothing,ftol=1.0/80,
    max_iter=10)

    @assert ndims(out_vol.vol)==3
    @assert ndims(in_vol.vol)==3

    if !isnothing(tfm) && isnothing(itfm)
        itfm=Minc2.inv(tfm)
    elseif isnothing(itfm) && isnothing(tfm)
        itfm=GeoTransforms() # identity transform
    end

    if isnothing(interp)
        if eltype(out_vol.vol) <: Integer # output is integer, use nearest neighbor
            interp=BSpline(Constant())
            if !isnothing(order)
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

    if isnothing(fill)
        fill=zero(eltype(out_vol.vol))
    end

    #TODO: extend this to grid support 
    resample_volume(in_vol.vol, out_vol.vol,in_vol.v2w, Minc2.inv(out_vol.v2w), itfm; interp,ftol, max_iter)
end

#TODO: make this more generic , and make it work with label 
function resample_grid_volume!(
    out_vol::Volume3D, in_vol::Volume3D; 
    tfm=nothing, itfm=nothing, interp=nothing, fill=nothing, 
    order=nothing,ftol=1.0/80,
    max_iter=10)

    @assert ndims(out_vol.vol)==4
    @assert ndims(in_vol.vol)==4

    if !isnothing(tfm) && isnothing(itfm)
        itfm=Minc2.inv(tfm)
    elseif isnothing(itfm) && isnothing(tfm)
        itfm=GeoTransforms() # identity transform
    end

    if isnothing(interp)
        if eltype(out_vol.vol) <: Integer # output is integer, use nearest neighbor
            interp=BSpline(Constant())
            if !isnothing(order)
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

    if isnothing(fill)
        fill=zero(eltype(out_vol.vol))
    end

    resample_grid_volume(in_vol.vol, out_vol.vol,in_vol.v2w, Minc2.inv(out_vol.v2w), itfm; interp,ftol, max_iter)
end
