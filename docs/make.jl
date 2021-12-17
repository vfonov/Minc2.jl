using Documenter
using Minc2

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "Minc2.jl Documentation",
    pages = [
        "Index" => "index.md",
        "Second page" => "secondPage.md"
    ]
    format = Documenter.HTML(prettyurls = false),
    modules = [Minc2]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
