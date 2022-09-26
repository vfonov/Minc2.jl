module Minc2

# only needed for world_to_voxel, matrix inversion
# and for factorization functions (?)
using LinearAlgebra

# geometric trasnsformation functions
include("geo_transform.jl")

# low level C-library interfaces to libminc
include("minc2_simple.jl")
# low level functions for interfacing with minc files
include("minc2_io.jl")
# low level functions for interfacing with .xfm files
include("xfm_io.jl")
# high level functions
include("minc_hl.jl")


end # module
