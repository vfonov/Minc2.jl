

# precompile some commonly used cases
precompile(read_ants_transform,(String,))

precompile(Volume3D,(Array{Float32,3}, AffineTransform{Float64},Nothing))
precompile(Volume3D,(Array{Float64,3}, AffineTransform{Float64},Nothing))
precompile(Volume3D,(Array{Int8,3},    AffineTransform{Float64},Nothing))
precompile(Volume3D,(Array{UInt8,3},    AffineTransform{Float64},Nothing))

precompile(empty_volume_like,(Volume3D{Float32,3},Float32,Nothing))
precompile(empty_volume_like,(Volume3D{Float64,3},Float64,Nothing))
precompile(empty_volume_like,(Volume3D{Int8,3},Int8,Nothing))
precompile(empty_volume_like,(Volume3D{UInt8,3},UInt8,Nothing))

precompile(save_volume,(String,Volume3D{Float32,3},Float32,Nothing))
precompile(save_volume,(String,Volume3D{Float32,3},Int16,Nothing))
precompile(save_volume,(String,Volume3D{Float64,3},Float32,Nothing))
precompile(save_volume,(String,Volume3D{Float64,3},Int16,Nothing))
precompile(save_volume,(String,Volume3D{Int8,3},Int8,Nothing))
precompile(save_volume,(String,Volume3D{UInt8,3},UInt8,Nothing))

precompile(read_volume,(String, Float64))
precompile(read_volume,(String, Float32))
precompile(read_volume,(String, Int8))
precompile(read_volume,(String, UInt8))

precompile(read_nifti_volume,(String,Float64))
precompile(read_nifti_volume,(String,Float32))
precompile(read_nifti_volume,(String,UInt8))

precompile(save_nifti_volume,(String,Volume3D{Float32,3},Float32,Nothing))
precompile(save_nifti_volume,(String,Volume3D{Float64,3},Float64,Nothing))
precompile(save_nifti_volume,(String,Volume3D{UInt8,3},UInt8,Nothing))

precompile(load_transforms,(String,))
precompile(save_transforms,(String,Vector{AnyTransform},Float32))
precompile(save_transforms,(String,Vector{AnyTransform},Float64))

precompile(transform_point,(IdentityTransform,SVector{3,Float64}))
precompile(transform_point,(AffineTransform{Float64},SVector{3,Float64}))
precompile(transform_point,(AffineTransform{Float64},SVector{3,Float64}))

#
GridTransform(Float64,Float64)
GridTransform(Float64,Float32)
InverseGridTransform(Float64,Float64)
InverseGridTransform(Float64,Float32)

