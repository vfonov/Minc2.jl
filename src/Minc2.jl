module Minc2

    using CBinding
    using LIBMINC_jll
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


end # module
