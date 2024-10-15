using CBinding
using .minc2_simple

"""
Low Level: Axis types from MINC volume
"""
@enum DIM begin
    DIM_UNKNOWN = Cint(minc2_simple.MINC2_DIM_UNKNOWN)
    DIM_X    = Cint(minc2_simple.MINC2_DIM_X )
    DIM_Y    = Cint(minc2_simple.MINC2_DIM_Y )
    DIM_Z    = Cint(minc2_simple.MINC2_DIM_Z )
    DIM_TIME = Cint(minc2_simple.MINC2_DIM_TIME )
    DIM_VEC  = Cint(minc2_simple.MINC2_DIM_VEC )
    DIM_END  = Cint(minc2_simple.MINC2_DIM_END )
end


"""
Low Level: Mapping MINC2 spatial dimensions to proper spatial dimension (should be identity for spatial dims)
"""
minc2_spatial=Dict(DIM_X=>1,
               DIM_Y=>2,
               DIM_Z=>3)

"""
Low Level: map Julia types to minc2 data types
"""
julia_to_minc2 = Dict(
    Type{Int8}     => Cint(minc2_simple.MINC2_BYTE ),
    Type{Int16}    => Cint(minc2_simple.MINC2_SHORT),
    Type{Int32}    => Cint(minc2_simple.MINC2_INT ),
    #Type{Int64}    => Cint(minc2_simple.MINC2_LONG ),
    Type{Float32}  => Cint(minc2_simple.MINC2_FLOAT ),
    Type{Float64}  => Cint(minc2_simple.MINC2_DOUBLE ),
    Type{String}   => Cint(minc2_simple.MINC2_STRING ),
    Type{UInt8}    => Cint(minc2_simple.MINC2_UBYTE ),
    Type{UInt16}   => Cint(minc2_simple.MINC2_USHORT ),
    Type{UInt32}   => Cint(minc2_simple.MINC2_UINT ),
    #Type{UInt64}   => Cint(minc2_simple.MINC2_ULONG ),
    Type{Complex{Int16}} => Cint(minc2_simple.MINC2_SCOMPLEX ),
    Type{Complex{Int32}} => Cint(minc2_simple.MINC2_ICOMPLEX ),
    Type{Complex{Float32}} => Cint(minc2_simple.MINC2_FCOMPLEX ),
    Type{Complex{Float64}} => Cint(minc2_simple.MINC2_DCOMPLEX ),
    Type{Any}      => Cint(minc2_simple.MINC2_UNKNOWN)
)

const _minc2_dimension=c"minc2_simple.struct minc2_dimension"


"""
Low Level: map MINC2 types to Julia
"""
minc2_to_julia=Dict([(j,i) for (i,j) in julia_to_minc2])

"""
Low Level: minc2_simple API status
"""
@enum STATUS begin
    # minc2 status
    SUCCESS  = Cint(minc2_simple.MINC2_SUCCESS)
    ERROR    = Cint(minc2_simple.MINC2_ERROR)
end

"""
MINC io exception
"""
struct Minc2Error <: Exception
    message::String
end

Base.showerror(io::IO, e::Minc2Error) = print(io,"MINC2:", e.message)

"""
    minc2_check( ex )

Macro to verify the return code
"""
macro minc2_check( ex ) # STATUS::SUCCESS
    return :($(esc(ex)) == 0 ? $(nothing) : throw( Minc2Error("error")) )
end

"""
    VolumeHandle

minc2_simple volume handle
"""
mutable struct VolumeHandle
    x::Ref
    function VolumeHandle()
        ret = new( Ref(minc2_simple.minc2_allocate0()) )

        finalizer(x -> c"minc2_simple.minc2_destroy"(x[]), ret.x)
        return ret
    end
end



"""
    MincHeader

Structure describing spatial orientation and sampling of the minc file 

"""
mutable struct MincHeader
    " Dimensions in voxels"
    dims::Vector{Int32}
    " Start offset (not origin)"
    start::Vector{Float64}
    " Sampling step"
    step::Vector{Float64}
    " Basis vector (direction cosine)"
    dir_cos::Matrix{Float64}
    " Flag that dir_cos contains valid information"
    dir_cos_valid::Vector{Bool}
    # TODO: figure out how to make it compatible with other libraries
    " Axis ID, see [DIM]"
    axis::Vector{DIM}

    """
    Construct header with given number of dimensions
    """
    function MincHeader(ndim)
        new(zeros(ndim),
            zeros(ndim), 
            ones(ndim), 
            Matrix{Float64}(undef, ndim, 3),
            zeros(Bool,ndim),
            fill(DIM_UNKNOWN,ndim))
    end
end

"""
    open_minc_file(fname::String)::VolumeHandle

Open minc file, return handle
"""
function open_minc_file(fname::String)::VolumeHandle
    h = VolumeHandle()
    @minc2_check minc2_simple.minc2_open(h.x[], fname)
    return h
end

"""
    define_minc_file(hdr::MincHeader,::Type{Store}=Float32,::Type{Repr}=Store) 

Define new minc file structure, return handle
"""
function define_minc_file(hdr::MincHeader,::Type{Store}=Float32,::Type{Repr}=Store) where {Store,Repr}
    h = VolumeHandle()

    _store_type = julia_to_minc2[Type{Store}]
    _representation_type = julia_to_minc2[Type{Repr}]

    _dims = [
        c"minc2_simple.struct minc2_dimension"(
                id=Int32(hdr.axis[i]),
                length=j,
                start=hdr.start[i],
                step=hdr.step[i],
                have_dir_cos=Int(hdr.dir_cos_valid[i]),
                dir_cos=Tuple(hdr.dir_cos[i,:])
            ) for (i,j) in enumerate(hdr.dims)
    ]
    push!(_dims,c"minc2_simple.struct minc2_dimension"(id=Cint(minc2_simple.MINC2_DIM_END)))

    @minc2_check minc2_simple.minc2_define(h.x[], _dims, _store_type, _representation_type)

    # TODO: fix this
    _global_scaling = 0
    _slice_scaling  = 0
    if _store_type != _representation_type && # HACK
            !(Store==Float32 && Repr==Float64 || Store==Float64 && Repr==Float32) &&
            !(Store==UInt8 && Repr==Int8 || Store==Int8 && Repr==UInt8) &&
            !(Store==UInt16 && Repr==Int16 || Store==Int16 && Repr==UInt16) 
        _global_scaling = 1 # slice scaling is not completely supported yet in minc2-simple
        #_slice_scaling = 1
    end
    @minc2_check minc2_simple.minc2_set_scaling(h.x[],_global_scaling,_slice_scaling )

    return h
end

"""
    create_minc_file(h::VolumeHandle, path::AbstractString)

Create empty minc file on disk, structure need to be defined in advance
"""
function create_minc_file(h::VolumeHandle, path::AbstractString)
    @minc2_check minc2_simple.minc2_create(h.x[], path )
end

"""
    close_minc_file(h::VolumeHandle)

Close currently open minc file, commit data on disk
"""
function close_minc_file(h::VolumeHandle)
    @minc2_check minc2_simple.minc2_close( h.x[] )
end

"""
    ndim(h::VolumeHandle)::Int

Query number of dimensions in the minc file
"""
function ndim(h::VolumeHandle)::Int
    dd = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_ndim( h.x[], dd )
    return dd[]
end

"""
    store_type(h::VolumeHandle)::Type

Query storage data type of the minc file, WARNING: it doesn't nacessery mean that it's the same as representation type
"""
function store_type(h::VolumeHandle)::Type
    minc_dtype = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_storage_data_type( h.x[], minc_dtype )
    return minc2_to_julia[minc_dtype[]]
end

"""
    store_type(h::VolumeHandle)::Type

Query representation data type of the minc file, WARNING: it doesn't nacessery mean that it's the same as store_type type
"""
function representation_type(h::VolumeHandle)::Type
    minc_dtype = Ref{Int}(0)
    @minc2_check minc2_simple.minc2_data_type( h.x[], minc_dtype )
    return minc2_to_julia[minc_dtype[]]
end


"""
    setup_standard_order(h::VolumeHandle)

Prepare to read volume in standard order: [V,X,Y,Z,TIME]
"""
function setup_standard_order(h::VolumeHandle)
    @minc2_check minc2_simple.minc2_setup_standard_order(h.x[])
end

"""
    _hdr_convert!(hdr::MincHeader,dd)::MincHeader

internal function
"""
function _hdr_convert!(hdr::MincHeader,dd)::MincHeader
    for i = 1:length(hdr.start)
        hdr.dims[i] = dd[][i].length
        hdr.start[i] = dd[][i].start
        hdr.step[i] = dd[][i].step
        if dd[][i].have_dir_cos != 0
            hdr.dir_cos[i,:] .= dd[][i].dir_cos[:]
            hdr.dir_cos_valid[i] = true
        else # TODO: do it only for spatial dimension
            hdr.dir_cos[i,:] = zeros(3)
            if haskey(minc2_spatial, DIM(dd[][i].id) )
                hdr.dir_cos[i, minc2_spatial[DIM(dd[][i].id)] ] = 1.0
            end
            hdr.dir_cos_valid[i] = false
        end
        hdr.axis[i] = DIM(dd[][i].id)
    end
    return hdr
end

"""
    representation_header(h::VolumeHandle)::MincHeader

Return volume representation header
"""
function representation_header(h::VolumeHandle)::MincHeader
    hdr = MincHeader(ndim(h))
    dd = Ref{c"minc2_simple.struct minc2_dimension *"}()

    # this is supposed to be in a standard order
    @minc2_check minc2_simple.minc2_get_representation_dimensions(h.x[], dd)

    return _hdr_convert!(hdr,dd)
end

"""
    store_header(h::VolumeHandle)::MincHeader

Return volume on-disk stucture header
"""
function store_header(h::VolumeHandle)::MincHeader
    hdr = MincHeader(ndim(h))
    dd = Ref{c"minc2_simple.struct minc2_dimension *"}()

    # this is supposed to be in a standard order
    @minc2_check minc2_simple.minc2_get_store_dimensions(h.x[], dd)

    return _hdr_convert!(hdr,dd)
end





"""
    empty_like_minc_volume_raw( h::VolumeHandle,
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader}

Allocate empty volume using handle
return volume, storage header
"""
function empty_like_minc_volume_raw( h::VolumeHandle,
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader} where {T<:Number}
    # TODO: use ImageMetadata to store header contents?
    store_hdr = store_header( h )
    volume = Array{T}(undef, store_hdr.dims...)

    return volume, store_hdr
end



"""
    Internal function that makes shure types are compatible
"""
function minc2_load_and_convert_complete_volume(
    h::VolumeHandle, volume::Array{T}, dtype::Type{T}) where {T}

    @minc2_check minc2_simple.minc2_load_complete_volume(h.x[], 
        Base.unsafe_convert(Ptr{Cvoid},volume), julia_to_minc2[Type{T}] )

    return volume
end

"""
    Internal function that makes sure types are compatible
"""
function minc2_load_and_convert_complete_volume(
    h::VolumeHandle, volume::Array{Tin}, dtype::Type{Tout}) where {Tin<:AbstractFloat,Tout<:Integer}

    @warn "Truncating volume data to integer type " Tin Tout

    @minc2_check minc2_simple.minc2_load_complete_volume(h.x[], 
        Base.unsafe_convert(Ptr{Cvoid},volume), julia_to_minc2[Type{Tin}] )

    return trunc.(Tout,volume)
end


"""
    Internal function that makes shure types are compatible
"""
function minc2_load_and_convert_complete_volume(
    h::VolumeHandle, volume::Array{Tin}, dtype::Type{Tout}) where {Tin<:Number,Tout<:Number}
    # make sure that types are compatible
    @minc2_check minc2_simple.minc2_load_complete_volume(h.x[], 
        Base.unsafe_convert(Ptr{Cvoid},volume), julia_to_minc2[Type{Tin}] )

    return convert.(Tout,volume)
end





"""
    read_minc_volume_raw(h::VolumeHandle, 
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader}
        
Read the actual volume using handle
return volume, storage header
"""
function read_minc_volume_raw(h::VolumeHandle,
        dtype::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader} where {T<:Number}

    fdtype=representation_type(h)
    volume, store_hdr = empty_like_minc_volume_raw(h, fdtype.parameters[1])

    return minc2_load_and_convert_complete_volume(h, volume, dtype), store_hdr
end


"""
    empty_like_minc_volume_std(h::VolumeHandle, 
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader}

Read the actual volume using handle
return volume, representation header,storage header
"""
function empty_like_minc_volume_std(h::VolumeHandle, 
    ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader} where {T<:Number}
    # TODO: use ImageMetadata to store header contents?
    setup_standard_order( h )
    store_hdr = store_header( h )
    hdr = representation_header( h )
    volume = Array{T}(undef, hdr.dims...)

    return volume, hdr, store_hdr
end


"""
    read_minc_volume_std(h::VolumeHandle, 
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader}
        
Read the actual volume using handle
return volume, representation header,storage header
"""
function read_minc_volume_std(h::VolumeHandle, 
        dtype::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader} where {T<:Number}
    # TODO: use ImageMetadata to store header contents?

    fdtype=representation_type(h)

    volume, hdr, store_hdr = empty_like_minc_volume_std(h,fdtype.parameters[1])

    return minc2_load_and_convert_complete_volume(h, volume, dtype), hdr,store_hdr

end


"""
    empty_like_minc_volume_std(path::String, 
        ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader}

allocate empty volume using path as a reference
return volume, representation header,storage header
"""
function empty_like_minc_volume_std(path::String, 
    ::Type{T}=Float32 )::Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader} where {T<:Number}
    handle = open_minc_file(path)
    volume, hdr, store_hdr = empty_like_minc_volume_std(handle,T)
    close_minc_file(handle)
    finalize(handle)
    return volume, hdr, store_hdr
end

"""
    empty_like_minc_volume_std_history(path::String, ::Type{T}=Float32 )

allocate empty volume using path
return volume, representation header,storage header
"""
function empty_like_minc_volume_std_history(path::String, ::Type{T}=Float32 ) where {T<:Number}
    handle = open_minc_file(path)
    volume, hdr, store_hdr = empty_like_minc_volume_std(handle,T)
    history = read_history(handle)
    close_minc_file(handle)
    finalize(handle)
    return volume, hdr, store_hdr, history
end


"""
    read_minc_volume_std(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader}

Read the actual volume using path
return volume, representation header,storage header
"""
function read_minc_volume_std(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader} where {T<:Number}
    handle = open_minc_file(path)
    volume, hdr, store_hdr = read_minc_volume_std(handle,T)
    close_minc_file(handle)
    finalize(handle)
    return volume, hdr, store_hdr
end


"""
    read_minc_volume_std_history(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader, Union{String,Nothing}}

Read the actual volume using path
return volume, representation header,storage header
"""
function read_minc_volume_std_history(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader, Minc2.MincHeader, Union{String,Nothing}} where {T<:Number} 
    handle = open_minc_file(path)
    volume, hdr, store_hdr = read_minc_volume_std(handle,T)
    history = read_history(handle)
    close_minc_file(handle)
    finalize(handle)
    return volume, hdr, store_hdr, history
end



"""
    empty_like_minc_volume_raw(path::String, ::Type{T}=Float32 )::
            Tuple{Array{T}, Minc2.MincHeader}

Create empty volume similar to existing file
return volume, representation header,storage header
"""
function empty_like_minc_volume_raw(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader} where {T<:Number}
    handle = open_minc_file(path)
    volume, store_hdr = empty_like_minc_volume_raw(handle,T)
    close_minc_file(handle)
    finalize(handle)
    return volume, store_hdr
end


"""
    read_minc_volume_raw(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader}

Read the actual volume using path
return volume, representation header,storage header
"""
function read_minc_volume_raw(path::String, ::Type{T}=Float32 )::
        Tuple{Array{T}, Minc2.MincHeader} where {T<:Number}
    handle = open_minc_file(path)
    volume, store_hdr = read_minc_volume_raw(handle,T)
    close_minc_file(handle)
    finalize(handle)
    return volume, store_hdr
end

"""
    read_minc_volume_raw_history(path::String, ::Type{T}=Float32 )

Read the actual volume using path
return volume, representation header,storage header
"""
function read_minc_volume_raw_history(path::String, ::Type{T}=Float32 ) where {T<:Number}
    handle = open_minc_file(path)
    volume, store_hdr = read_minc_volume_raw(handle,T)
    history = read_history(h)
    close_minc_file(handle)
    finalize(handle)
    return volume, store_hdr, history
end



"""
    write_minc_volume_raw(h::VolumeHandle, volume::Array{T} )

write full volume to file, file should be defined and created
return nothing
"""
function write_minc_volume_raw(h::VolumeHandle, volume::Array{T} ) where {T<:Number}
    @minc2_check minc2_simple.minc2_save_complete_volume(h.x[], Base.unsafe_convert(Ptr{Cvoid},volume),julia_to_minc2[Type{T}])
    return nothing
end


"""
    write_minc_volume_std(h::VolumeHandle, volume::Array{T} ) 

write full volume to file, file should be defined and created
return nothing
"""
function write_minc_volume_std(h::VolumeHandle, volume::Array{T} ) where {T<:Number}
    setup_standard_order( h )
    write_minc_volume_raw(h, volume)
    return nothing
end

"""
    copy_minc_metadata(i::VolumeHandle, o::VolumeHandle)

Copy metadata from one file to another
"""
function copy_minc_metadata(i::VolumeHandle, o::VolumeHandle)
    @minc2_check minc2_simple.minc2_copy_metadata(i.x[],o.x[])
end


"""
    read_attribute(h::VolumeHandle, 
        group::String,
        attribute::String)::Union{String, AbstractVector, Nothing}

Read minc2 header attribute
:param group: attribute group name
:param attribute: attribute name
:return:
"""
function read_attribute(h::VolumeHandle, 
    group::String,
    attribute::String)::Union{String, AbstractVector, Nothing}

    attr_type  = Ref{Int}(0)
    attr_length = Ref{Int}(0)

    # assume that if we can't get attribute type, it's missing, return nil

    @minc2_check minc2_simple.minc2_get_attribute_type(h.x[], group, attribute, attr_type)
    @minc2_check minc2_simple.minc2_get_attribute_length(h.x[], group, attribute, attr_length)

    # assume that array of BYTES is a string 
    # TODO: make sure it's a good assumption
    if attr_type[] == minc2_simple.MINC2_STRING 
        # remove trailing 0
        _buf = Vector{UInt8}(undef,attr_length[]-1) 
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, _buf, attr_length[])
        return String(_buf)
    elseif attr_type[] == minc2_simple.MINC2_FLOAT
        # TODO: finish this
        buf = Vector{Float32}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_DOUBLE
        # TODO: finish this
        buf = Vector{Float64}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_INT
        # TODO: finish this
        buf = Vector{Int32}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_SHORT
        # TODO: finish this
        buf = Vector{Int16}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_UINT
        # TODO: finish this
        buf = Vector{UInt32}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_USHORT
        # TODO: finish this
        buf = Vector{UInt16}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_BYTE
        # TODO: finish this
        buf = Vector{Int8}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    elseif attr_type[] == minc2_simple.MINC2_UBYTE
        # TODO: finish this
        buf = Vector{UInt8}(undef,attr_length[])
        @minc2_check minc2_simple.minc2_read_attribute(h.x[], group, attribute, buf, attr_length[])
        return buf
    else
        throw(SystemError("MINC2 unsupported data type"))
    end

    return nothing
end

"""
    get_attribute(h::VolumeHandle,g::String,a::String; default=missing)

Convenience function for reading specific attribute, return default value if not found
also convert Array into the first value if it's a one-length array
"""
function get_attribute(h::VolumeHandle,g::String,a::String; default=missing)
    # ::Union{String,Missing,...}
    r=default
    if g in groups(h)
        if a in attributes(h,g)
            r=read_attribute(h,g,a)
        end
    end

    if !ismissing(r)
        if typeof(r) <: Array{T, 1} where T<:Number
            if length(r)==1
                r=r[1]
            end
        end
        if typeof(r) <: String
            r=rstrip(r,'\0') # HACK
        end
    end
    return r
end


"""
    groups(h::VolumeHandle)

List groups defined in minc2 file
"""
function groups(h::VolumeHandle)
    i = Ref(minc2_simple.minc2_allocate_info_iterator())
    r = Vector{String}()
    try
        @minc2_check minc2_simple.minc2_start_group_iterator(h.x[],i[])
        while minc2_simple.minc2_iterator_group_next(i[]) == minc2_simple.MINC2_SUCCESS
            push!(r,unsafe_string(minc2_simple.minc2_iterator_group_name(i[])))
        end
    finally
        @minc2_check minc2_simple.minc2_free_info_iterator(i[])
    end
    return r
end


"""
    attributes(h::VolumeHandle, g::String)

List attributes defined inside given group
"""
function attributes(h::VolumeHandle, g::String)
    i = Ref(minc2_simple.minc2_allocate_info_iterator())
    r = Vector{String}()
    try
        @minc2_check minc2_simple.minc2_start_attribute_iterator(h.x[],g,i[])
        while minc2_simple.minc2_iterator_attribute_next(i[]) == minc2_simple.MINC2_SUCCESS
            push!(r,unsafe_string(minc2_simple.minc2_iterator_attribute_name(i[])))
        end
    finally
        @minc2_check minc2_simple.minc2_free_info_iterator(i[])
    end
    return r
end

"""
    write_attribute(h::VolumeHandle, group::String, attribute::String, value::T)

Store attribute into minc2 file
:param group:  group name
:param attribute:  attribute name
:param value:  attribute value
:return:
"""
function write_attribute(h::VolumeHandle, group::String, attribute::String, value::T) where {T}
    val_type = julia_to_minc2[Type{T}]
    attr_length = length(value)
    # TODO: make this work for other type
    # Right now it might write garbage for anything but strings
    if isa(value,String)
        _value=Vector{UInt8}(value)
        @minc2_check minc2_simple.minc2_write_attribute(h.x[], group, attribute, Base.unsafe_convert(Ptr{Cvoid},_value), attr_length, val_type)
    else
        @minc2_check minc2_simple.minc2_write_attribute(h.x[], group, attribute, Base.unsafe_convert(Ptr{Cvoid},value), attr_length, val_type)
    end
end

"""
    read_history(i::VolumeHandle)::Union{String, Nothing}

Return history string
"""
function read_history(i::VolumeHandle)::Union{String, Nothing}
    return get_attribute(i, "", "history", default=nothing)
end

"""
    write_history(i::VolumeHandle, history::String)

Write history string
"""
function write_history(i::VolumeHandle, history::String)
    write_attribute(i,"","history",history)
end

"""
    world_to_voxel(h::VolumeHandle, xyz::Vector{Float64})::Vector{Float64}

Convert world coordinates (X,Y,Z) to contignuous voxel indexes (I,J,K) 0-based
"""
function world_to_voxel(h::VolumeHandle, xyz::Vector{Float64})::Vector{Float64}
    # TODO: check input vector size (?)
    @assert(length(xyz)==3)
    ijk::Vector{Float64}=zeros(Float64,3)
    @minc2_check minc2_simple.minc2_world_to_voxel(h.x[],xyz,ijk)
    return ijk
end

"""
    voxel_to_world(hdr::MincHeader)::AffineTransform{Float64}

Give AffineTransform for world to voxel transformation based on header
"""
function voxel_to_world(hdr::MincHeader)::AffineTransform{Float64}
    rot=zeros(3,3)
    scales=zeros(3,3)
    start=zeros(3)
    for i=1:3
        aa=findfirst(isequal(DIM(i)), hdr.axis) # HACK, assumes DIM_X=1,DIM_Y=2 etc
        if hdr.dir_cos_valid[aa]
            rot[1:3,i] = hdr.dir_cos[aa,1:3]
        else
            rot[i,i] = 1
        end
        scales[i,i] = hdr.step[aa]
        start[i] = hdr.start[aa]
    end
    origin = rot*start 

    return AffineTransform(rot*scales, origin)
end

"""
    world_to_voxel(hdr::MincHeader)::AffineTransform{Float64}

Give AffineTransform for voxel to world transformation
"""
function world_to_voxel(hdr::MincHeader)::AffineTransform{Float64}
    inv(voxel_to_world(hdr))
end

"""
    create_header_from_v2w(
        sz, t::AffineTransform{T};
        vector_dim::Bool=false, 
        time_step::Union{Float64,Nothing}=nothing,
        time_start::Union{Float64,Nothing}=nothing)::MincHeader

Internal: Generate header from the voxel to world transform and volume size
"""
function create_header_from_v2w(
        sz, t::AffineTransform{T};
        vector_dim::Bool=false, 
        time_step::Union{Float64,Nothing}=nothing,
        time_start::Union{Float64,Nothing}=nothing)::MincHeader where {T}

    time_dim = !isnothing(time_step) && !isnothing(time_start)

    start, step, dir_cos = decompose(t)

    hdr = MincHeader(3 + vector_dim + time_dim )

    # TODO: verify sz dimensionality
    for i = 1:length(sz)
        hdr.dims[i] = sz[i]

        if vector_dim && i==1 # vector dimension
            hdr.start[i] = 0
            hdr.step[i] = 1.0
            hdr.dir_cos_valid[i] = false
            hdr.dir_cos[i,:] = [0.0,0.0,0.0]
            hdr.axis[i] = DIM_VEC
        elseif time_dim && (i-vector_dim)==4 # time dimension
            hdr.step[i] = time_step
            hdr.start[i] = time_start
            hdr.axis[i] = DIM_TIME
        else # spatial dimension
            hdr.start[i] = start[i-vector_dim]
            hdr.step[i] = step[i-vector_dim]
            hdr.dir_cos[i,:] .= dir_cos[:, i-vector_dim]
            hdr.dir_cos_valid[i] = true
            hdr.axis[i] = DIM(i-vector_dim)
        end
    end
    return hdr
end



"""
    voxel_to_world(h::VolumeHandle,ijk::Vector{Float64})::Vector{Float64}

Convert contignuous 0-based voxel indexes (I,J,K) to world coordinates (X,Y,Z) 0-based
"""
function voxel_to_world(h::VolumeHandle,ijk::Vector{Float64})::Vector{Float64}
    @assert(length(ijk)==3)
    xyz::Vector{Float64}=zeros(Float64,3)
    @minc2_check minc2_simple.minc2_voxel_to_world(h.x[],ijk,xyz)
    return xyz
end


"""
    write_minc_volume_std(path::String, ::Type{Store}, 
        store_hdr::Union{MincHeader, Nothing}, 
        volume::Array{Repr};like::Union{String, Nothing}=nothing, 
        history::Union{String, Nothing}=nothing )

Write full volume to file, need to provide details of file structure
return nothing
"""
function write_minc_volume_std(path::String, ::Type{Store}, 
    store_hdr::Union{MincHeader, Nothing}, 
    volume::Array{Repr};like::Union{String, Nothing}=nothing, 
    history::Union{String, Nothing}=nothing ) where {Store,Repr}

    if isnothing(like)
        # TODO: check if store_hdr is compatible with volume

        handle = define_minc_file(store_hdr, Store, Repr)
        create_minc_file(handle,path)

        if !isnothing(history)
            write_history(handle,history)
        end

        write_minc_volume_std(handle,volume)
        close_minc_file(handle)
    else
        # need to copy file structure
        in_h = open_minc_file(like)
        store_hdr = store_header( in_h )

        # TODO: check if store_hdr is compatible with volume
        handle = define_minc_file(store_hdr, Store, Repr)
        create_minc_file(handle,path)
        
        copy_minc_metadata(in_h,handle)

        if !isnothing(history)
            write_history(handle,history)
        end

        write_minc_volume_std(handle,volume)
        close_minc_file(in_h)
        close_minc_file(handle)
    end
    finalize(handle)

    return nothing
end

"""
    write_minc_volume_raw(path::String, ::Type{Store}, 
        store_hdr::Union{MincHeader, Nothing}, 
        volume::Array{Repr};like::Union{String, Nothing}=nothing, 
        history::Union{String, Nothing}=nothing )

Write full volume to file, need to provide details of file structure
return nothing
"""
function write_minc_volume_raw(path::String, ::Type{Store}, 
    store_hdr::Union{MincHeader, Nothing}, 
    volume::Array{Repr};like::Union{String, Nothing}=nothing, 
    history::Union{String, Nothing}=nothing ) where {Store, Repr}

    if isnothing(like)
        # TODO: check if store_hdr is compatible with volume

        handle = define_minc_file(store_hdr, Store, Repr)
        create_minc_file(handle,path)

        if !isnothing(history)
            write_minc_history(handle,history)
        end

        write_minc_volume_raw(handle, volume)
        close_minc_file(handle)
    else
        # need to copy file structure
        in_h = open_minc_file(like)
        store_hdr = store_header( in_h )

        # TODO: check if store_hdr is compatible with volume
        handle = define_minc_file(store_hdr, Store, Repr)
        create_minc_file(handle,path)
        
        copy_minc_metadata(in_h,handle)

        if !isnothing(history)
            write_minc_history(handle,history)
        end

        write_minc_volume_raw(handle,volume)
        close_minc_file(in_h)
        close_minc_file(handle)
    end
    finalize(handle)
    return nothing
end
