# Filippov Example

This is an example example. These examples could probably somehow be automatically populated & updated from the Pluto example file.

```@example
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings

F(x) = [3, -1]
G(x) = [0, 1]
H(x) = x[2] - sin(x[1])
N(x) = [-cos(x[1]), 1]

sysF = HD.FilippovSystem(F, G, H, N)
probF = HD.prob(sysF, [0.0, 1.0], (0.0, 10.0))
solF = HD.solve(probF, HD.RK4())

xf = getindex.(solF.x, 1)
yf = getindex.(solF.x, 2)
xh = range(minimum(xf)-0.5, maximum(xf)+0.5, length=1_000)

plt.plot(xf, yf, lw=2, label="Filippov Trajectory")
plt.plot!(xh, sin.(xh), lw=2, lc=:black, ls=:dash, label=L"H(x)=0")
plt.plot!(title = "Filippov Trajectory", label=L"x", ylabel=L"y")
```