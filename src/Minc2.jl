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

    struct Volume
        x::Ref
        function Volume()
          ret = new(Ref(minc2_allocate0()))
          finalizer(x -> minc2_destroy(x[]), ret.x)
          return ret
        end
    end
    


end # module
