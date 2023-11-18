"""
Collection of methods for reading and writing MINC2 files,
    and for reading and writing .xfm files.
    and operations on 3D volumes
"""
module Minc2

# only needed for world_to_voxel, matrix inversion
# and for factorization functions (?)
using LinearAlgebra
using Requires

# geometric transformation functions
include("geo_transform.jl")
# low level C-library interfaces to libminc
include("minc2_simple.jl")
# low level functions for interfacing with minc files
include("minc2_io.jl")
# low level functions for interfacing with .xfm files
include("xfm_io.jl")
# high level functions
include("minc_hl.jl")
# ITK and NIFTI interface
include("nifti_io.jl")

# draw illustration using Makie, if available
function __init__()
    @require Colors="5ae59095-9a9b-59fe-a467-6f913c188581" begin
        @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
            import .Colors
            import .Makie

            include("minc_makie.jl")

            # Drawing using Makie
            export draw_outline_with_labels
        end
    end
end


# Try to precompile things
include("precompile.jl")

end # module
