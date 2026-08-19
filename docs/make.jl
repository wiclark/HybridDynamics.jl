using Documenter
using HybridDynamics

# Automatically populate index.md with the git readme
write(joinpath(@__DIR__, "src", "index.md"),
      read(joinpath(".", "README.md"), String))

# Make the build
makedocs(
    sitename = "HybridDynamics",
    format = Documenter.HTML(), 
    modules = [HybridDynamics],
    checkdocs = :none,              # Remove when...documentation is written?
    # remotes = nothing,               # Remove when...the git repo is public?
    pages = [
        "Home" => "index.md",
        "System Types" => Any[
            "General" => "System_Types/General.md",
            "Linear and Affine" => "System_Types/Linear_Affine.md",
            "Mechanical" => "System_Types/Mechanical.md",
            "Nonholonomic" => "System_Types/Nonholonomic.md",
            "Stochastic" => "System_Types/Stochastic.md",
            "Filippov" => "System_Types/Filippov.md"
        ],
        "Analysis" => Any[
            "Utility Functions" => "Analysis/Utilities.md",
            "Lyapunov Exponents" => "Analysis/Lyapunov_exp.md"
        ],
        "Solver Algorithms" => "Solver_Algorithms/all.md",
        "Examples" => Any[
            "Basic Systems" => Any[
                "General" => "Examples/Basic_sys/ex_general.md",
                "Linear and Affine" => "Examples/Basic_sys/ex_linear_affine.md",
                "Mechanical" => "Examples/Basic_sys/ex_mechanical.md",
                "Nonholonomic" => "Examples/Basic_sys/ex_nonholonomic.md",
                "Stochastic" => "Examples/Basic_sys/ex_stochastic.md",
                "Filippov" => "Examples/Basic_sys/ex_filippov.md"
            ],
            "Advanced" => Any[
                "Lyapunov Exponents" => "Examples/Advanced/ex_lyapunovexp.md"
            ]
        ],
        "Change Log" => "changelog.md",
        "Known Bugs" => "knownissues.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/wiclark/HybridDynamics.jl.git",
    devbranch = "main"
)
