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

    struct VolumeHandle
        x::Ref
        function VolumeHandle()
          ret = new(Ref(c"minc2_simple.minc2_allocate0"()))
          finalizer(x -> c"minc2_simple.minc2_destroy"(x[]), ret.x)
          return ret
        end
    end

    mutable struct minc_header
        start::Vector{Float64,1}
        step::Vector{Float64}
        dir_cos::Array{Float64,2}
        # history?
        # tags?
    end

    function open_minc_file(fname::String)
        h::VolumeHandle
        c"minc2_simple.minc2_open"(fname,h.x)
        return h
    end

    function close_minc_file(h::VolumeHandle)
        c"minc2_simple.minc2_close"(h.x)
    end

    function setup_standard_order(h::VolumeHandle)
        c"minc2_simple.setup_standard_order"(h.x)
    end
    
    function read_minc_volume(h::VolumeHandle)
        dims=c"minc2_simple.minc2_dimension"()
        lib.minc2_get_representation_dimensions(h.x,dims)
    end


end # module
