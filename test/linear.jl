# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 0.5

e = 0.8
    
A = [0.0 1.0 0.0;
    0.0 0.0 -9.81;
    0.0 0.0 0.0]

λ = [1.0, 0.0, 0.0]

C = [1.0 0.0 0.0;
    0.0 -e 0.0;
    0.0 0.0 1.0]

x0 = [10.0, 0.0, 1.0]

sysL = HD.LinearSystem(A, λ, C)
probL = HD.prob(sysL, x0, (0.0, 2.0))

@testset "Linear" begin
    
    solL = HD.solve(probL, HD.RK4())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solL(1) ≈ [5.095, -9.81, 1] atol=soltol
    # After bounce, position
    @test solL(1.75)[1] ≈ -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)) atol=soltol
    # One bounce
    @test length(solL.event_times) == 1
    # At the right time
    @test solL.event_times[1] ≈ sqrt(10*2/9.81) atol=soltol

end