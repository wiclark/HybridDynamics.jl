# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 0.5

# Define problem with known solution
# This is the (holonomic) pendulum disguised as a nonholonomic system.
MNH(q) = [1.0 0.0; 0.0 1.0]
VNH(q) = q[2]
ANH(q) = [q[1] q[2]]
hNH(q) = q[1]

sysNH = HD.NonholonomicSystem(MNH, VNH; A=ANH, guard=hNH, e=1.0)
probNH = HD.prob(sysNH, [1.0, 0.0, 0.0, 0.0], (0.0, 4.0))

@testset "Nonholonomic" begin
    
    solNH = HD.solve(probNH, HD.RK4())

    # Test to verify that a single bounced occurred
    @test length(solNH.event_times) == 1

    # Test to ensure the constraint is satisfied
    @test ANH(solNH.x[end][1:2])*solNH.x[end][3:4] ≈ [0.0] atol=soltol
end