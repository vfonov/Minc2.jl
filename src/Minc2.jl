module Minc2

# only needed for world_to_voxel, matrix inversion
using LinearAlgebra

include("geo_tools.jl")

# low leverl C-library interfaces
include("minc2_simple.jl")
include("minc2_io.jl")
include("xfm_io.jl")


end # module
