using Clang.Generators
using Pkg.Artifacts
using Base.BinaryPlatforms

cd(@__DIR__)

# Find libminc_jll's Artifacts.toml via pkgdir
# (import instead of using — avoids triggering artifact_dir resolution)
import libminc_jll
artifact_toml = joinpath(pkgdir(libminc_jll), "Artifacts.toml")
@assert isfile(artifact_toml) "Could not find Artifacts.toml at $artifact_toml"

# Build platform with MPI tag
platform = HostPlatform()
platform["mpi"] = "mpitrampoline"

h = artifact_hash("libminc", artifact_toml; platform=platform)
@assert h !== nothing "Could not find artifact hash for libminc on platform $platform"

include_dir = joinpath(artifact_path(h), "include")
@assert isdir(include_dir) "Include directory not found: $include_dir"

@info "Using include directory: $include_dir"

# Generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# Compiler args (include paths)
args = get_default_args()
push!(args, "-I$include_dir")

# The header file to parse
headers = [joinpath(include_dir, "minc2-simple.h")]
@assert isfile(headers[1]) "Header not found: $(headers[1])"

# Create and run generator
@info "Generating bindings for minc2-simple..."
ctx = create_context(headers, args, options)
build!(ctx)

@info "Done!"
