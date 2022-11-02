# for resampling
using Interpolations
# for minc history
using Dates

"""
An abstract 3D volume
"""
struct Volume3D
    vol # an abstract volume
    v2w # transformation from voxel to world coordinates
    history # file metadata: history
end

function Volume3D(vol, v2w::AffineTransform{T}; history=nothing) where {T}
    return Volume3D(vol, v2w, history)
end


function Volume3D(vol, like::Volume3D; history=nothing)
    return Volume3D(vol, like.v2w, isnothing(history) ? like.history : history)
end

function read_volume(fn::String; store::Type{T}=Float64) where {T}
    in_vol,in_hdr,in_store_hdr,in_history = Minc2.read_minc_volume_std_history(fn, store)
    v2w=Minc2.voxel_to_world(in_hdr)

    return Volume3D(in_vol,v2w,in_history)
end


function empty_volume_like(fn::String; store::Type{T}=Float64, history=nothing) where {T}
    out_vol,out_hdr, out_store_hdr, ref_history = Minc2.empty_like_minc_volume_std_history(fn,store)
    v2w=Minc2.voxel_to_world(out_hdr)

    return Volume3D(out_vol,v2w,isnothing(history) ? ref_history : history)
end


function empty_volume_like(vol::Volume3D; store::Type{T}=Float64, history=nothing) where {T}
    out_vol = Array{Float64}(undef, size(vol.vol)...)
    return Volume3D(out_vol, vol.v2w, isnothing(history) ? vol.history : history )
end


function save_volume(fn, vol::Volume3D; store::Type{T}=Float32,history=nothing) where {T}
    if isnothing(history)
        _history=vol.history
    else
        _history=(isnothing(vol.history) ? "" : vol.history * "\n" )*history
    end

    Minc2.write_minc_volume_std(fn, store, 
      Minc2.create_header_from_v2w(size(vol.vol), vol.v2w, 
      vector_dim=(length(size(vol.vol))==4)), vol.vol; history=_history)
end



function resample_volume!(in_vol::Array{T,3}, 
    out_vol::Array{T,3}, 
    v2w::Minc2.AffineTransform{C}, 
    w2v::Minc2.AffineTransform{C}, 
    itfm::Vector{XFM};
    interp::I=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10) where {C, T, I, XFM<:Minc2.AnyTransform}

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
function resample_grid_volume!(
    in_vol::Array{T,4},
    out_vol::Array{T,4},
    v2w::Minc2.AffineTransform{C}, 
    w2v::Minc2.AffineTransform{C}, 
    itfm::Vector{XFM};
    interp::I=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10) where {C, T, I, XFM<:Minc2.AnyTransform}

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
function resample_grid(in_grid, itfm; like=nothing)::Minc2.Volume3D
    if isnothing(like)
      out_vol = similar(in_grid.vol)
      v2w = in_grid.v2w
    else
      out_vol = similar(like.vol)
      v2w = like.v2w
    end
    resample_grid_volume!(in_grid.vol, out_vol, in_grid.v2w, Minc2.inv(v2w), itfm; interp=BSpline(Linear()))
    return Minc2.Volume3D(out_vol, v2w)
end

function tfm_to_grid!(
        tfm::Vector{XFM}, 
        grid::G,
        v2w::Minc2.AffineTransform{C}) where {T,C, XFM<:Minc2.AnyTransform, G<:AbstractArray}
    
    @simd for c in CartesianIndices(size(grid)[2:end])
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(tfm, orig)
        @inbounds grid[:,c] .= dst .- orig
    end
end

# convert transforms into a single nonlinear grid
function tfm_to_grid(tfm::Vector{XFM},
        ref::G; 
        store::Type{T}=Float64)::Minc2.Volume3D where {T, XFM<:Minc2.AnyTransform, G<:Minc2.GridTransform}
    out_grid = similar(ref.vector_field, store)
    v2w = ref.voxel_to_world

    tfm_to_grid!(tfm,out_grid,v2w)
    return Minc2.Volume3D( out_grid, v2w)
end

# convert transforms into a single nonlinear grid transform
function normalize_tfm(tfm::Vector{XFM},
    ref::G; 
    store::Type{T}=Float64)::Minc2.GridTransform{Float64,T} where {T, XFM<:Minc2.AnyTransform, G<:Minc2.GridTransform}

    out_grid = similar(ref.vector_field, store)
    v2w = ref.voxel_to_world

    @simd for c in CartesianIndices(size(out_grid)[2:end])
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(tfm, orig)

        @inbounds out_grid[:,c] .= dst .- orig
    end

    return Minc2.GridTransform{Float64,T}(v2w, out_grid)
end

#TODO: make this more generic , and make it work with labels 
function resample_volume!(in_vol::Volume3D, out_vol::Volume3D; 
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
    resample_volume!(in_vol.vol, out_vol.vol, out_vol.v2w, Minc2.inv(in_vol.v2w), itfm; 
        interp,ftol, max_iter)

    return out_vol
end


function resample_volume(in_vol::Volume3D;like=nothing,
    tfm=nothing, itfm=nothing, interp=nothing, fill=nothing, 
    order=nothing,ftol=1.0/80, store=nothing,
    max_iter=10)

    out_vol=Minc2.empty_volume_like(isnothing(like) ? in_vol : like ;
        store = isnothing(store) ? (isnothing(like) ? eltype(in_vol.vol) : eltype(like.vol) ) : store)

    return resample_volume!(in_vol,out_vol;tfm,itfm,interp,fill,order,ftol,max_iter)
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

    resample_grid_volume!(in_vol.vol, out_vol.vol,in_vol.v2w, Minc2.inv(out_vol.v2w), itfm; interp,ftol, max_iter)
end

function calculate_jacobian!(
    out_vol::Array{T,3}, 
    v2w::Minc2.AffineTransform{C},
    tfm::GeoTransforms;
    interp::I=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10) where {C,T,I}

    # calculate scaling matrix from the voxel to world matrix
    _,step,_ = decompose(v2w)
    sc = diagm(step) #Base.inv(v2w.rot * Base.inv(dir_cos))
    #@info "Scaling matrix" sc

    # First step: generate vector field of transformations
    vector_field = Array{T}(undef, 3, size(out_vol)...)

    @simd for c in CartesianIndices(out_vol)
        orig = Minc2.transform_point(v2w, c )
        dst  = Minc2.transform_point(tfm, orig; ftol, max_iter )

        @inbounds vector_field[:,c] .= dst # .- orig
    end

    # Second step: calculate jacobian determinant 
    vector_field_itp = extrapolate( interpolate( vector_field, (NoInterp(),interp,interp,interp)), Flat())

    @simd for c in CartesianIndices(out_vol)
        grad = hcat([ Interpolations.gradient( vector_field_itp, i, Tuple(c)...) for i in 1:3 ]...)
        out_vol[c] = det(grad'*sc)
    end

    out_vol
end

function calculate_jacobian!(out_vol::Volume3D, tfm::GeoTransforms; interp=BSpline(Quadratic(Line(OnCell()))),
    fill=0.0,
    ftol=1.0/80,
    max_iter=10)
    calculate_jacobian!(out_vol.vol,out_vol.v2w,tfm;interp,fill,ftol,max_iter)
end

# generate minc-style history from program args
function format_history(args)
    #stamp=strftime("%a %b %d %T %Y>>>", gmtime())
    # Thu Jul 30 14:23:47 2009
    return Dates.format(now(),"e d u HH:MM:SS YYYY")*">>>"*join(args," ")
end
