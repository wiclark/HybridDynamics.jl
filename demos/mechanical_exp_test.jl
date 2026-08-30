import HybridDynamics as HD
import ForwardDiff
using LinearAlgebra

# A ball in a disk
M(q) = [1.0 0.0; 0.0 1.0]
V(q) = q[2]
h(q) = 1 - (q[1]^2+q[2]^2)
∇h(q) = [-2*q[1], -2*q[2]]
r = 1.0

sysM = HD.MechanicalSystem(M, V; guard=h, normal=∇h, e=r)
probM = HD.prob(sysM, [0.2, 0.0, 0.0, 0.0], (0.0, 2.0))
solM = HD.solve(probM, HD.RK4(); tol=1e-9, dt_initial=1e-3)

# The STM and exponents
ΦM = HD.tangent_dynamics(solM, sysM)
λ = HD.LyapunovExponents(sysM, [0.2, 0.0, 0.0, 0.0]; run_length=1.0, run_iter=Int(1e3), transient=1_000.0, dt_initial=1e-2, solver=HD.RK4())