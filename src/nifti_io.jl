using NIfTI


function read_nifti_volume(fn::String; store::Type{T}=Float64) where {T}
    ni = niread(fn)

    #convert volume to a simple volume and afine matrix to the v2w 
    # check if the scaling is making sense
    if !(ni.header.scl_slope==0.0f0 && ni.header.scl_inter==0.0f0)
        in_vol=map(x->convert(store, x),ni)
    else # we can't trust scaling
        in_vol=map(x->convert(store, x),ni.raw)
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

    start, step, dir_cos = decompose(vol.v2w)
    orientation=Matrix{Float32}([dir_cos start'])

    ni=NIVolume(map(x->convert(store, x),vol.vol); 
                voxel_size=Tuple(Float32(x) for x in step),
                orientation= orientation)
    # TODO: deal with vectors (?)
    niwrite(fn,ni)
end


function read_ants_warp(fn::String)
    V=read_nifti_volume(fn; store=Float64)
    # need to reshape to confirm to the MINC convention 
    return V # TODO: finish this
end

function read_itk_transform(fn::String)
    # replicating https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/IO/TransformInsightLegacy/src/itkTxtTransformIO.cxx
    transform_type=""
    parameters=Float64[]
    fixed_parameters=Float64[]
    transform_type=""
    for line in eachline(fn)
        line = strip(line,[' ', '\t', '\n','\r'])

        if length(line) == 0
              continue
        elseif line[1]=='#'
              # Skip lines beginning with #, or blank lines
              continue
        end
        # parse tags: name : value 
        # will throw an error if format is different
        @info line
        (tag_name, tag_val) = split(line,":"; limit=2)
        tag_name = strip(tag_name,[' ', '\t', '\n','\r'])
        if tag_name == "Transform"
            transform_type = strip(tag_val,[' ', '\t', '\n','\r'])
        elseif tag_name == "ComponentTransformFile"
            # not supported for now
        elseif tag_name == "Parameters"
            # read transform paramters
            parameters = parse.(Float64,split(tag_val," ",keepempty=false))
        elseif tag_name == "FixedParameters"
            # read transform fixed parameters
            fixed_parameters = parse.(Float64,split(tag_val," ",keepempty=false))
        else
            # something we don't know about, ignore ?
            @warn "Uknown tag" tag_name
        end
    end

    # deal with transformations ? 
    (tfm_type, tfm_data_type, tfm_dims) = split(transform_type,'_'; limit=3)
    if tfm_type == "MatrixOffsetTransformBase" || tfm_type == "AffineTransform"
        # ignoring fixed parameters
        @assert tfm_dims=="3_3"
        @assert length(parameters) == 12
        return Minc2.AffineTransform(reshape(parameters[1:9],3,3),parameters[10:12])
    else
        throw( Minc2Error("Unsupported transform type \"$tfm_type\""))
    end
end
