# List of Solver Algorithms

### Linear Multi-step Methods

#### Fixed:
```@docs
AdamsBashforth2
AdamsBashforth3
BDF2
```

#### Adaptive:
```@docs
AdaptiveABM2
AdaptiveABM3
```

---
### Runge-Kutta Methods

#### Fixed:
```@autodocs
Modules = [HybridDynamics]
Pages = ["ODE_Solvers/Runge_Kutta_Methods/Fixed_RK_Steps.jl"]
```

#### Adaptive:
```@autodocs
Modules = [HybridDynamics]
Pages = ["ODE_Solvers/Runge_Kutta_Methods/Adaptive_RK_Steps.jl"]
```

---
### Miscelaneous
```@docs
ExponentialSolver
EulerMaruyama
MagnusLeapfrog
```