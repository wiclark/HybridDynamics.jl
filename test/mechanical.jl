# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 0.5

# Define problem with known solution
M(q) = 1
V(q) = q[1]*9.81
h(q) = q[1]

sysM = HD.MechanicalSystem(M, V; guard=h, e=0.8)
probM = HD.prob(sysM, [10.0, 0.0], (0.0, 2.0))

@testset "Mechanical" begin
    
    solM = HD.solve(probM, HD.RK4())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solM(1) ≈ [5.095, -9.81] atol=soltol
    # After bounce, position
    @test solM(1.75)[1] ≈ -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)) atol=soltol
    # One bounce
    @test length(solM.event_times) == 1
    # At the right time
    @test solM.event_times[1] ≈ sqrt(10*2/9.81) atol=soltol

end