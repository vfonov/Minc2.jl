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

