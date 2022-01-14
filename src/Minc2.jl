module Minc2

    module minc2_simple
        using LIBMINC_jll
        using CBinding
        # Main.SYSROOT...,
        c`$([
            "-I$(LIBMINC_jll.artifact_dir)/include",
            "-L$(dirname(LIBMINC_jll.libminc2_path))", 
            "-L$(dirname(LIBMINC_jll.libminc2_simple_path))", 
            "-llibminc2-simple","-llibminc2",
        ])`

        c"""
        #include <minc2-simple.h>
        """ji

        #c"typedef double minc2_dir_cos[3];"j
    end

    using CBinding
    using .minc2_simple
    using LinearAlgebra


    @enum DIM begin
        DIM_UNKNOWN = Cint(minc2_simple.MINC2_DIM_UNKNOWN)
        DIM_X = Cint(minc2_simple.MINC2_DIM_X  )
        DIM_Y = Cint(minc2_simple.MINC2_DIM_Y )
        DIM_Z = Cint(minc2_simple.MINC2_DIM_Z )
        DIM_TIME = Cint(minc2_simple.MINC2_DIM_TIME )
        DIM_VEC  = Cint(minc2_simple.MINC2_DIM_VEC )
        DIM_END  = Cint(minc2_simple.MINC2_DIM_END)
    end

    # map Julia types to minc2 data types
    julia_to_minc2 = Dict(
       Type{Int8} => Cint(minc2_simple.MINC2_BYTE ),
       Type{Int16}    =>Cint(minc2_simple.MINC2_SHORT),
       Type{Int32}      => Cint(minc2_simple.MINC2_INT ),
       Type{Float32}    => Cint(minc2_simple.MINC2_FLOAT ),
       Type{Float64}   => Cint(minc2_simple.MINC2_DOUBLE ),
       Type{String}   => Cint(minc2_simple.MINC2_STRING ),
       Type{UInt8}    => Cint(minc2_simple.MINC2_UBYTE ),
       Type{UInt16}   => Cint(minc2_simple.MINC2_USHORT ),
       Type{UInt32}     => Cint(minc2_simple.MINC2_UINT ),
       Type{Complex{Int16}} => Cint(minc2_simple.MINC2_SCOMPLEX ),
       Type{Complex{Int32}} => Cint(minc2_simple.MINC2_ICOMPLEX ),
       Type{Complex{Float32}} => Cint(minc2_simple.MINC2_FCOMPLEX ),
       Type{Complex{Float64}} => Cint(minc2_simple.MINC2_DCOMPLEX ),
       Type{Any}  => Cint(minc2_simple.MINC2_UNKNOWN)
    )

    const _minc2_dimension=c"minc2_simple.struct minc2_dimension"
    

    # map MINC2 types to Julia
    minc2_to_julia=Dict([(j,i) for (i,j) in julia_to_minc2])

    @enum STATUS begin
        # minc2 status
        SUCCESS  = Cint(minc2_simple.MINC2_SUCCESS)
        ERROR    = Cint(minc2_simple.MINC2_ERROR)
    end

    macro minc2_check( ex ) # STATUS::SUCCESS
        return :($(esc(ex)) == 0 ? $(nothing) : throw(SystemError("MINC2 error")))
    end

    mutable struct VolumeHandle
        x::Ref
        function VolumeHandle()
          ret = new( Ref(minc2_simple.minc2_allocate0()) )

          finalizer(x -> c"minc2_simple.minc2_destroy"(x[]), ret.x)
          return ret
        end
    end

    mutable struct MincHeader
        dims::Vector{Int32}
        start::Vector{Float64}
        step::Vector{Float64}
        dir_cos::Matrix{Float64}

        # TODO: figure out how to make it compatible with other libraries
        axis::Vector{Int32}

        function MincHeader(ndim)
            new(zeros(ndim), zeros(ndim), ones(ndim), Matrix{Float64}(I, 3, 3), zeros(ndim))
        end
    end

    function open_minc_file(fname::String)
        h = VolumeHandle()
        @minc2_check minc2_simple.minc2_open(h.x[], fname)
        return h
    end

    function define_minc_file(hdr::MincHeader,::Type{A}=Float32,::Type{B}=A) where {A,B}
        h = VolumeHandle()

        _store_type = julia_to_minc2[Type{A}]
        _representation_type = julia_to_minc2[Type{B}]
        
        _dims = [
            c"minc2_simple.struct minc2_dimension"(
                    id=hdr.axis[i],
                    length=j,
                    start=hdr.start[i],
                    step=hdr.step[i],
                    have_dir_cos=0,
                    dir_cos=Tuple(hdr.dir_cos)
                ) for (i,j) in enumerate(hdr.dims)
        ]
        push!(_dims,c"minc2_simple.struct minc2_dimension"(id=Cint(minc2_simple.MINC2_DIM_END)))

        @minc2_check minc2_simple.minc2_define(h.x[], _dims, _store_type, _representation_type)

        # TODO: fix this
        _global_scaling = 0
        _slice_scaling  = 0 
        if _store_type != _representation_type
            _global_scaling = 1 # slice scaling is not completely supported yet in minc2-simple
            #_slice_scaling = 1
        end
        @minc2_check minc2_simple.minc2_set_scaling(h.x[],_global_scaling,_slice_scaling )

        return h
    end

    function create_minc_file(h::VolumeHandle, path::String)
        @minc2_check minc2_simple.minc2_create(h.x[], path )
    end

    function close_minc_file(h::VolumeHandle)
        @minc2_check minc2_simple.minc2_close( h.x[] )
    end

    function ndim(h::VolumeHandle)
        dd = Ref{Int}(0)
        @minc2_check minc2_simple.minc2_ndim( h.x[], dd )
        return dd[]
    end

    function setup_standard_order(h::VolumeHandle)
        @minc2_check minc2_simple.minc2_setup_standard_order(h.x[])
    end
    
    function representation_header(h::VolumeHandle)::MincHeader
        hdr = MincHeader(ndim(h))
        dd = Ref{c"minc2_simple.struct minc2_dimension *"}()

        # this is supposed to be in a standard order
        @minc2_check minc2_simple.minc2_get_representation_dimensions(h.x[], dd)

        for i = 1:length(hdr.start)
            hdr.dims[i] = dd[][i].length
            hdr.start[i] = dd[][i].start
            hdr.step[i] = dd[][i].step
            if dd[][i].have_dir_cos != 0
                hdr.dir_cos[i,:] .= dd[][i].dir_cos[:]
            hdr.axis[i] = dd[][i].id
            end
        end
        
        return hdr
    end

    function store_header(h::VolumeHandle)::MincHeader
        hdr = MincHeader(ndim(h))
        dd = Ref{c"minc2_simple.struct minc2_dimension *"}()

        # this is supposed to be in a standard order
        @minc2_check minc2_simple.minc2_get_store_dimensions(h.x[], dd)

        for i = 1:length(hdr.start)
            hdr.dims[i] = dd[][i].length
            hdr.start[i] = dd[][i].start
            hdr.step[i] = dd[][i].step
            if dd[][i].have_dir_cos != 0
                hdr.dir_cos[i,:] .= dd[][i].dir_cos[:]
            end
            hdr.axis[i] = dd[][i].id
        end
        
        return hdr
    end

    function read_minc_volume(h::VolumeHandle, ::Type{T}=Float32 ) where {T}
        # TODO: use ImageMetadata to store header contents?
        setup_standard_order( h )
        store_hdr = store_header( h )
        hdr = representation_header( h )
        volume = Array{T}(undef, hdr.dims...)

        @minc2_check minc2_simple.minc2_load_complete_volume(h.x[], Base.unsafe_convert(Ptr{Cvoid},volume), julia_to_minc2[Type{T}] )
        return volume, hdr, store_hdr
    end

    function write_minc_volume(h::VolumeHandle, volume::Array{T} ) where {T}
        setup_standard_order( h )
        @minc2_check minc2_simple.minc2_save_complete_volume(h.x[], Base.unsafe_convert(Ptr{Cvoid},volume),julia_to_minc2[Type{T}])
    end


end # module
