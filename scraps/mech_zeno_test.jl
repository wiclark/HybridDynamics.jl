import HybridDynamics as HD
import Plots as plt

# I want to see how our predicted Zeno time compares to the exact for the bouncing ball.

# The theoretical value
zeno_time(x₀, p₀, e, g) = (p₀+sqrt(2*g*x₀+sqrt(p₀^2+2*g*x₀)))/g + 2*e*sqrt(p₀^2+2*g*x₀)/(g*(1-e))

g, e = 9.81, 0.5
M(q) = [1.0;;]
V(q) = g*q[1]
h(q) = q[1]
∇h(q) = [1.0]

sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=e)
probM = HD.prob(sysM, [10.0, 0.0], (0.0, 15.0))
solM = HD.solve(probM, solver=HD.RK45(); ztol=1e-5)

tz = solM.zeno[1]
println(tz)
println(zeno_time(10.0, 0.0, e, g))

## Let's do some testing
function zeno_tolerence(z_tol)
    solM = HD.solve(probM, solver=HD.RK4(); ztol=z_tol, dt_initial=1e-4)
    return abs(zeno_time(10.0, 0.0, e, g) - solM.zeno[1])
end

Ztol  = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
EZeno = zeno_tolerence.(Ztol)
plt.scatter(Ztol, EZeno, xaxis=:log, label="")