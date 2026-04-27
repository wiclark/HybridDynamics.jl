# I want to understand hybrid Lagrangian systems (Hamiltonian will follow)

struct HybridLagrangianSystem
    L::Function
    h::Function
    e::Function
end

# I want the following:
# 1. Generate trajectories (you can call forward diff)
# 2. Implement an intelligent way to perform event detection
# 3. Determine whether or not the system is Zeno and state when/where that occurs.