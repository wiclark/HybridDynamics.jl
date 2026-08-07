# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 0.5

g = 9.81
e = 0.8

A = [0.0 1.0;
     0.0 0.0]

b = [0.0, -g]

λ = [1.0, 0.0]

a = 0.0

C = [1.0 0.0;
     0.0 -e]

κ = [0.0, 0.0]

x0 = [10.0, 0.0]

sysA = HD.AffineSystem(A, b, λ, a, C, κ)
probA = HD.prob(sysA, x0, (0.0, 2.0))

@testset "Affine" begin
    
    solA = HD.solve(probA, HD.RK4())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solA(1) ≈ [5.095, -9.81] atol=soltol
    # After bounce, position
    @test solA(1.75)[1] ≈ -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)) atol=soltol
    # One bounce
    @test length(solA.event_times) == 1
    # At the right time
    @test solA.event_times[1] ≈ sqrt(10*2/9.81) atol=soltol

end