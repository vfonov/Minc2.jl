# interface with minc2-simple C library
# linked with libminc 
module minc2_simple

using libminc_jll
using CBinding

c`$([
    "-I$(libminc_jll.artifact_dir)/include",
    "-L$(dirname(libminc_jll.libminc2_path))", 
    "-L$(dirname(libminc_jll.libminc2_simple_path))", 
    "-llibminc2-simple","-llibminc2",
])`

c"""
#include <minc2-simple.h>
"""ji

end
