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

    @enum TYPE begin
        #: minc2 data types
       BYTE     = Cint(minc2_simple.MINC2_BYTE )
       SHORT    = Cint(minc2_simple.MINC2_SHORT)
       INT      = Cint(minc2_simple.MINC2_INT )
       FLOAT    = Cint(minc2_simple.MINC2_FLOAT )
       DOUBLE   = Cint(minc2_simple.MINC2_DOUBLE )
       STRING   = Cint(minc2_simple.MINC2_STRING )
       UBYTE    = Cint(minc2_simple.MINC2_UBYTE )
       USHORT   = Cint(minc2_simple.MINC2_USHORT )
       UINT     = Cint(minc2_simple.MINC2_UINT )
       SCOMPLEX = Cint(minc2_simple.MINC2_SCOMPLEX )
       ICOMPLEX = Cint(minc2_simple.MINC2_ICOMPLEX )
       FCOMPLEX = Cint(minc2_simple.MINC2_FCOMPLEX )
       DCOMPLEX = Cint(minc2_simple.MINC2_DCOMPLEX )
       MAX_TYPE_ID = Cint(minc2_simple.MINC2_MAX_TYPE_ID)
       UNKNOWN  = Cint(minc2_simple.MINC2_UNKNOWN)
    end

    @enum STATUS begin
        # minc2 status
        SUCCESS  = Cint(minc2_simple.MINC2_SUCCESS)
        ERROR    = Cint(minc2_simple.MINC2_ERROR)
    end

    macro minc2_check( ex )
        return :($(esc(ex)) == 0 ? $(nothing) : throw(SystemError("MINC2 error in " * __module__ * " " *__source__)))
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

        function MincHeader(ndim)
            new(zeros(ndim),zeros(ndim),ones(ndim),Matrix{Float64}(I, 3, 3))
        end
    end

    function open_minc_file(fname::String)
        h = VolumeHandle()
        @minc2_check minc2_simple.minc2_open(h.x[], fname)
        return h
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
    
    function header(h::VolumeHandle)
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
            end
        end
        
        return hdr
    end


    function read_minc_volume(h::VolumeHandle)
        setup_standard_order( h )
        hdr = header( h )
        volume = Array{Float32}(undef,hdr.dims...)
        @minc2_check minc2_simple.minc2_load_complete_volume(h.x[], Base.unsafe_convert(Ptr{Cvoid},volume), minc2_simple.MINC2_FLOAT )
        return hdr,volume
    end


end # module
