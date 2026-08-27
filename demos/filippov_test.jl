import HybridDynamics as HD
import Plots as plt

F(x) = [3, -1]
G(x) = [0, 1]
H(x) = x[2] - sin(x[1])
N(x) = [-cos(x[1]), 1]

sysF = HD.FilippovSystem(F, G, H, N)
probF = HD.prob(sysF, [0.0, 1.0], (0.0, 10.0))
solF = HD.solve(probF, HD.RK45())

# Try out Lyapunov!
λ = HD.LyapunovExponents(sysF, [0.0, 1.0]; run_iter = 100, transient=100.0)