using Clang.Generators
using libminc_jll
using Pkg

cd(@__DIR__)

# Get include directory from jll
include_dir = joinpath(libminc_jll.artifact_dir, "include")

# Generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# Compiler args (include paths)
args = get_default_args()
push!(args, "-I$include_dir")

# The header file to parse
headers = [joinpath(include_dir, "minc2-simple.h")]

# Create and run generator
@info "Generating bindings for minc2-simple..."
ctx = create_context(headers, args, options)

# Write output to ../src/LibMinc2Simple.jl
build!(ctx)

