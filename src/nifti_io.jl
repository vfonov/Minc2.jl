using NIfTI
using Rotations
using StaticArrays

function read_nifti_volume(fn::String; store::Type{T}=Float64) where {T}
    ni = niread(fn)

    #convert volume to a simple volume and afine matrix to the v2w 
    # check if the scaling is making sense
    if !(ni.header.scl_slope==0.0f0 && ni.header.scl_inter==0.0f0)
        in_vol=map(x->convert(store, x),ni)
    else # we can't trust scaling
        in_vol=map(x->convert(store, x),ni.raw)
    end

    tfm = Float64.(getaffine(ni))

    # need to flip x,y 
    # to simulate behaviour of ITK
    # see https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/IO/NIFTI/src/itkNiftiImageIO.cxx#L2031
    tfm[1:3,1:2] .= tfm[1:3,1:2] .* -1.0
    tfm[1:2,4]   .= tfm[1:2,4]   .* -1.0

    v2w=Minc2.AffineTransform( tfm )

    return Minc2.Volume3D(in_vol, v2w, fn)
end


function save_nifti_volume(fn, vol::Volume3D; store::Type{T}=Float32,history=nothing) where {T}
    if isnothing(history)
        _history=vol.history
    else
        _history=vol.history*"\n"*history
    end
    tfm = [vol.v2w.rot vol.v2w.shift; 0 0 0 1]
    # flip 
    tfm[1:3,1:2] .= tfm[1:3,1:2] .* -1.0
    tfm[1:2,4]   .= tfm[1:2,4]   .* -1.0

    start, step, dir_cos = decompose(tfm)
    R_quat = Rotations.params(QuatRotation(dir_cos))

    ni=NIVolume(map(x->convert(store, x),vol.vol),
        qfac=1.0f0,
        quatern_b=  R_quat[2],  quatern_c=R_quat[3],  quatern_d=R_quat[4],
        qoffset_x= start[1,1], qoffset_y=start[1,2], qoffset_z=start[1,3],
        voxel_size=Tuple(Float32(i) for i in step),
        xyzt_units=Int8(2),
        regular=Int8('r'))
    
    #setaffine(ni.header, [vol.v2w.rot vol.v2w.shift;0 0 0 1])
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
    _delimiters=[' ', '\t', '\n','\r']
    parameters=Float64[]
    fixed_parameters=Float64[]
    transform_type=""
    for line in eachline(fn)
        line = strip(line,_delimiters)

        if length(line) == 0
              continue
        elseif line[1]=='#'
              # Skip lines beginning with #, or blank lines
              continue
        end
        # parse tags: name : value 
        # will throw an error if format is different
        #@info line
        (tag_name, tag_val) = split(line,":"; limit=2)
        tag_name = strip(tag_name,_delimiters)
        if tag_name == "Transform"
            transform_type = strip(tag_val,_delimiters)
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
