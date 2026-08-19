import HybridDynamics as HD
import ForwardDiff
import Plots as plt
using LinearAlgebra

#=
## The elastic ball
function f_ball(x, t)
	g = 9.81
	α = 0.0
	q, v = x
	return [v, -g-α*v*abs(v)]
end

h_ball(x) = x[1]
Δ_ball(x) = [abs(x[1]), -1*x[2]]

sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
probG = HD.prob(sysG, [10.0, 0.0], (0.0, 10.0))
solG = HD.solve(probG, HD.RK4(), dt_initial=1e-3)

# λG = LyapunovExponents(sysG, [10.0, 0.0])
=#
## Bouncing ball in a circle
function f_circle(z, t)
    x, y, px, py = z
    return [px, py, 0, -1]
end
h_circle(z) = 1 - (z[1]^2 + z[2]^2)
function Δ_circle(z)
    x, y, px, py = z
    pxp = px - (2*x*(px*x + py*y))/(x^2 + y^2)
    pyp = py - (2*y*(px*x + py*y))/(x^2 + y^2)
    return [x, y, pxp, pyp]
end

sysC = HD.GeneralSystem(f_circle, h_circle, Δ_circle; direction=-1)
probC = HD.prob(sysC, [-0.95, 0.0, 0.2, -1.5], (0.0, 10.0))
solC = HD.solve(probC, HD.RK4(), dt_initial=1e-3)

xC, yC = getindex.(solC.x, 1), getindex.(solC.x, 2)
plt.plot(xC, yC)
θ = LinRange(0, 2π, 100)
plt.plot!(cos.(θ), sin.(θ))

λC = HD.LyapunovExponents(sysC, [-0.95, 0.0, 0.2, -1.5]; run_length=5.0)
println(λC)

function exponent_circle(x, y)
    if x^2+y^2 >= 1
        return NaN
    else
        return maximum(HD.LyapunovExponents(sysC, [x, y, 0.0, 0.0]; run_length=10.0, run_iter=Int(1e3)))
    end
end

res = 4
X = LinRange(-1,1,res)
Y = LinRange(-1,1,res)

# Threads.nthreads()

Z = zeros(res, res)

@time Threads.@threads for i ∈ 1:res
    for j ∈ 1:res
        Z[i,j] = exponent_circle(X[i], Y[j])
    end
end
##
plt.heatmap(X, Y, Z', dpi = 200)

## The Lorenz system - a benchmark
#=
function f_lorenz(q, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    x, y, z = q
    return [σ*(y-x), x*(ρ-z)-y, x*y-β*z]
end
h_lorenz(x) = 10.0
Δ_lorenz(x) = []

sysL = HD.GeneralSystem(f_lorenz, h_lorenz, Δ_lorenz)
probL = HD.prob(sysL, [1.0, 0.0, 0.0], (0.0, 10.0))
solL = HD.solve(probL, HD.RK4(), dt_initial=1e-3)

# λL = LyapunovExponents(sysL, [1.0, 0.0, 0.0]; transient=200.0)
# println(λL)

x_sol, y_sol, z_sol = getindex.(solL.x, 1), getindex.(solL.x, 2), getindex.(solL.x, 3)
plt.plot(x_sol, y_sol, z_sol)
=#