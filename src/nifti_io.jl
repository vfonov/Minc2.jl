using NIfTI
using Rotations
using StaticArrays


"""
Convert coordinates from LPS to RAS
"""
function lps_to_ras(xfm::Minc2.AffineTransform)
    # convert from LPS to RAS
    # see
    lps_to_ras_m = SMatrix{3,3,Float64}([-1 0 0; 0 -1 0; 0 0 1])
    out=Minc2.AffineTransform(lps_to_ras_m*xfm.rot, lps_to_ras_m*xfm.shift)
end


"""
Convert coordinates from LPS to RAS
"""
function lps_to_ras(vol::Minc2.Volume3D{T})::Minc2.Volume3D{T} where {T}
    Minc2.Volume3D(vol.vol, lps_to_ras(vol.v2w), vol.history)
end


"""
Read Volume3D from .nii or .nii.gz file
"""
function read_nifti_volume(fn::AbstractString; store::Type{T}=Float64)::Volume3D{T} where {T}
    ni = niread(fn)

    #convert volume to a simple volume and afine matrix to the v2w 
    tfm = Float64.(getaffine(ni.header))

    # need to flip x,y 
    # to simulate behaviour of ITK
    # see https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/IO/NIFTI/src/itkNiftiImageIO.cxx#L2031
    tfm[1:3,1:2] .= tfm[1:3,1:2] .* -1.0
    tfm[1:2,4]   .= tfm[1:2,4]   .* -1.0

    v2w = Minc2.AffineTransform( tfm )

    return Minc2.Volume3D( ! ((ni.header.scl_slope == 0.0f0 && ni.header.scl_inter == 0.0f0) ||
                              (ni.header.scl_slope == 1.0f0 && ni.header.scl_inter == 0.0f0) ) ? 
                               convert(AbstractArray{T}, ni) : convert(AbstractArray{T}, ni.raw) , 
                          v2w, fn)
end


"""
Save Volume3D into .nii or .nii.gz file
"""
function save_nifti_volume(fn::AbstractString, vol::Volume3D{T}; 
        store::Type{S}=Float32, history=nothing) where {T,S}
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

    ni=NIVolume( convert(AbstractArray{S}, vol.vol),
        qfac=1.0f0,
        quatern_b=  R_quat[2],  quatern_c=R_quat[3],  quatern_d=R_quat[4],
        qoffset_x= start[1], qoffset_y=start[2], qoffset_z=start[3],
        voxel_size=Tuple(Float32(i) for i in step),
        xyzt_units=Int8(2),
        regular=Int8('r'))
    
    #setaffine(ni.header, [vol.v2w.rot vol.v2w.shift;0 0 0 1])
    # TODO: deal with vectors (?)
    niwrite(fn,ni)
end


"""
Read ANTs style warp transform
"""
function read_itk_nifti_transform(fn::AbstractString;
        store::Type{T}=Float32)::Minc2.GridTransform{Float64,T} where {T}
    V = read_nifti_volume(fn; store)
    
    @assert ndims(V.vol)  == 5
    @assert size(V.vol,5) == 3
    #reshape into minc convention
    return Minc2.GridTransform(V.v2w, permutedims(dropdims(V.vol, dims=4),(4,1,2,3)))
end


"""
Write ANTs style warp transform
"""
function write_itk_nifti_transform(fn::AbstractString,
        xfm::Minc2.GridTransform{Float64,T}; store::Type{S}=Float32) where {T,S}
    save_nifti_volume(fn,Volume3D( permutedims(xfm.vector_field[:,:,:,:,:],(2,3,4,5,1)), xfm.voxel_to_world);store)
end


"""
Read .txt and .nii(.nii.gz) transforms produces by ANTs
"""
function read_ants_transform(fn::AbstractString; store::Type{T}=Float32)::Minc2.AnyTransform where {T}
    if endswith(fn,".txt")
        return read_itk_txt_transform(fn)
    elseif endswith(fn,".nii") || endswith(fn,".nii.gz")
        return read_itk_nifti_transform(fn;store)
    else
        @warn "Not sure about this transform: " fn
        return read_itk_nifti_transform(fn;store)
    end
end


"""
Read ITK legacy transform in .txt format
"""
function read_itk_txt_transform(fn::AbstractString)::Minc2.AffineTransform{Float64}
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
            parameters = parse.(Float64, split(tag_val," ",keepempty=false))
        elseif tag_name == "FixedParameters"
            # read transform fixed parameters
            fixed_parameters = parse.(Float64, split(tag_val," ",keepempty=false))
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
        @assert length(fixed_parameters) == 3
        rot_matrix   = reshape(parameters[1:9],3,3)
        translation  = parameters[10:12]
        center       = fixed_parameters

        # need to convert to rot_matrix and offset
        # as defined in 
        # https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Core/Transform/include/itkMatrixOffsetTransformBase.hxx#L552
        offset = center + translation - (center' * rot_matrix)'

        return Minc2.AffineTransform(rot_matrix, offset)
    else
        throw( Minc2Error("Unsupported transform type \"$tfm_type\""))
    end
end


"""
Write affine transform in ITK .txt format
"""
function write_itk_txt_transform(fn::AbstractString,
        xfm::Minc2.AffineTransform{Float64})
    # TODO: preserve fixed parameters?
    fixed_parameters = zeros(3) 
    parameters = vcat(reshape(xfm.rot,9), xfm.shift)
    open(fn,"w") do io
        write(io,"#Insight Transform File V1.0\n")
        write(io,"#Transform 0\n")
        write(io,"Transform: MatrixOffsetTransformBase_double_3_3\n")
        write(io,"Parameters: ",join(parameters," "),"\n")
        write(io,"FixedParameters: ",join(fixed_parameters," "),"\n")
    end
end
