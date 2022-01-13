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


    # @enum DIM begin 
    #     DIM_UNKNOWN=minc2_simple.MINC2_DIM_UNKNOWN 
    #     DIM_X = minc2_simple.MINC2_DIM_X  
    #     DIM_Y = minc2_simple.MINC2_DIM_Y 
    #     DIM_Z = minc2_simple.MINC2_DIM_Z 
    #     DIM_TIME = minc2_simple.MINC2_DIM_TIME 
    #     DIM_VEC  = minc2_simple.MINC2_DIM_VEC 
    #     DIM_END  = minc2_simple.MINC2_DIM_END
    # end

    # @enum TYPE begin
    #     #: minc2 data types
    #    BYTE     = minc2_simple.MINC2_BYTE 
    #    SHORT    = minc2_simple.MINC2_SHORT
    #    INT      = minc2_simple.MINC2_INT 
    #    FLOAT    = minc2_simple.MINC2_FLOAT 
    #    DOUBLE   = minc2_simple.MINC2_DOUBLE 
    #    STRING   = minc2_simple.MINC2_STRING 
    #    UBYTE    = minc2_simple.MINC2_UBYTE 
    #    USHORT   = minc2_simple.MINC2_USHORT 
    #    UINT     = minc2_simple.MINC2_UINT 
    #    SCOMPLEX = minc2_simple.MINC2_SCOMPLEX 
    #    ICOMPLEX = minc2_simple.MINC2_ICOMPLEX 
    #    FCOMPLEX = minc2_simple.MINC2_FCOMPLEX 
    #    DCOMPLEX = minc2_simple.MINC2_DCOMPLEX 
    #    MAX_TYPE_ID = minc2_simple.MINC2_MAX_TYPE_ID
    #    UNKNOWN  = minc2_simple.UNKNOWN
    # end

    # @enum STATUS begin
    #     # minc2 status
    #     SUCCESS  = minc2_simple.MINC2_SUCCESS,
    #     ERROR    = minc2_simple.MINC2_ERROR
    # end

    macro minc2_check( ex )
        return :($(esc(ex)) == 0 ? $(nothing) : throw(SystemError("MINC2 error in ")))
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
        start::Vector{Float64}
        step::Vector{Float64}
        dir_cos::Matrix{Float64}

        function MincHeader(ndim)
            new(zeros(ndim),ones(ndim),Matrix{Float64}(I, 3, 3))
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
        dd=Ref{Int}(0)
        @minc2_check minc2_simple.minc2_ndim( h.x[], dd )
        return dd[]
    end

    function setup_standard_order(h::VolumeHandle)
        @minc2_check minc2_simple.setup_standard_order(h.x[])
    end
    
    function header(h::VolumeHandle)
        hdr = MincHeader(ndim(h))
        dd = Ref{c"minc2_simple.struct minc2_dimension *"}()

        @minc2_check minc2_simple.minc2_get_representation_dimensions(h.x[], dd)
        for i=1:length(hdr.start)
            hdr.start[i] = dd[][i].start
            hdr.step[i] = dd[][i].step
            if dd[][i].have_dir_cos != 0
                hdr.dir_cos[i,:] .= dd[][i].dir_cos[:]
            end
        end
        
        return hdr
    end

    function read_minc_volume(h::VolumeHandle)
        dims=c"minc2_simple.minc2_dimension"()
        hdr=header(h)

        @minc2_check minc2_simple.minc2_get_representation_dimensions(h.x[], dims)
        ndims::Integer=ndim(h)

        return ndims
    end


end # module
