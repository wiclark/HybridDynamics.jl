
struct GeneralSystem <: AbstractHybridSystem
    n::Int          #Dimension of system
    f::Function     #Continuous Dynamics: (x,t) -> dx/dt
    h::Function     #Guard Surface: x-> real
    Δ::Function     #Reset map: x-> x⁺
end

struct GeneralProblem <: AbstractHybridProblem
    sys::GeneralSystem
    x₀::Vector{Float64}
    tspan::Tuple{Float64, Float64}
end

struct GeneralSolution 
    t::Vector{Float64}
    x::Vector{Vector{Float64}}
    jump_times::Vector{Float64}
    jump_indices::Vector{Int}
end

get_dimension(sys::GeneralSystem) = sys.n
vector_field(sys::GeneralSystem) = sys.f
guard(sys::GeneralSystem, x) = sys.h(x)
apply_reset(sys::GeneralSystem, x) = sys.Δ(x)


function GeneralProblem(f::Function, h::Function, Δ::Function, x₀::Vector{Float64}, tspan)
    n = length(x₀)
    sys = GeneralSystem(n, f, h, Δ)
    return GeneralProblem(sys, x₀, tspan)
end

function GeneralProblem(f::Function, x₀::Vector{Float64}, tspan)
    n = length(x₀)
    dummy_h(x) = 1.0
    dummy_Δ(x) = x
    sys = GeneralSystem(n, f, dummy_h, dummy_Δ)
    return GeneralProblem(sys, x₀, tspan)
end

function init_solution(prob::GeneralProblem)
    return GeneralSolution([prob.tspan[1]], [prob.x₀], Float64[], Int[])
end
