using Documenter
using Minc2

makedocs(
    sitename = "Minc2.jl",
    format   = Documenter.HTML(size_threshold = nothing,),
    modules  = [Minc2],
    doctest = false,
    warnonly = true,
    clean    = true,
    pages    = [
        "Introduction to Minc2.jl" => "index.md",

        "Reference" => [
            "Alphabetical function list" => "reference/index.md",
            "Function reference"         => "reference/api.md"
        ],
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(;
     repo = "github.com/vfonov/Minc2.jl.git",
     versions = nothing # temporary
)
