using Documenter
using HybridDynamics

makedocs(
    sitename = "HybridDynamics",
    format = Documenter.HTML(),
    modules = [HybridDynamics]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
