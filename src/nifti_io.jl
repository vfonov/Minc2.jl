using NIfTI


function read_nifti_volume(fn::String; store::Type{T}=Float64)
    ni = niread(fn)

    #convert volume to a simple volume and afine matrix to the v2w 
    # check if the scaling is making sense
    if !(ni.header.scl_slope==0.0f0 && ni.header.scl_inter==0.0f0)
        in_vol=map(x->convert(store,x),ni)
    else # we can't trust scaling
        in_vol=map(x->convert(store,x.raw),ni)
    end
    v2w=Minc2.AffineTransform(Float64.(getaffine(ni)))

    return Minc2.Volume3D(in_vol, v2w, fn)
end



function save_nifti_volume(fn, vol::Volume3D; store::Type{T}=Float32,history=nothing) where {T}
    if isnothing(history)
        _history=vol.history
    else
        _history=vol.history*"\n"*history
    end

    
    ni=NIVolume(vol.vol; voxel_size=voxel_size,
    orientation=orientation, dim_info=dim_info,
    time_step=time_step != false && !isempty(time_step.data) ? time_step.data[1] : 0f0)
end
