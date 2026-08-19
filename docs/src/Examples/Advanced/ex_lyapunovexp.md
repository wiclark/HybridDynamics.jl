# Lyapunov Exponents Example
!!! info "Bouncing ball in a vertical circle"
    Consider the general hybrid system with continuous dynamics
    ```math
        \ddot{x} = 0, \ddot{y} = -1, \quad x^2+y^2 < 1 \\
    ```
    with reset
    ```math
        \begin{split}
            \dot{x}^+ = \dot{x} - 2\frac{x(x\dot{x}+y\dot{y})}{x^2+y^2} \\
            \dot{y}^+ = \dot{y} - 2\frac{y(x\dot{x}+y\dot{y})}{x^2+y^2}
        \end{split}
    ```
    The goal is to numerically compute the Lyapunov exponents for this system.

```@example lyap_exp
import HybridDynamics as HD
```

```@example lyap_exp
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

sys = HD.GeneralSystem(f_circle, h_circle, Δ_circle; direction=-1)
```

```@example lyap_exp
λ = HD.LyapunovExponents(sys, [-0.5,0.5,0.0,0.0]; run_length=5.0)
```