using Documenter
using Minc2

makedocs(
    sitename = "Minc2",
    format = Documenter.HTML(),
    modules = [Minc2],
    doctest = false

)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
