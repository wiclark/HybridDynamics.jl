# I want a hybrid linear system and hybrid affine structures
struct HybridLinearSystem
    A::Matrix
    λ::Vector
    C::Matrix
end

# I want the following:
# 1. Generate trajectories
# 2. Determine whether or not the system is trivially blocking
# 3. Find the beating/blocking sets
# 4. (For affine) determine whether or not the system is Zeno and find the Zeno time