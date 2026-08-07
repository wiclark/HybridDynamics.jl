# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 1e-2

# Define problem with known solution
function f_ball(x, t)
    g = 9.81
    q, v = x
    return [v, -g]
end
h_ball(x) = x[1]
Δ_ball(x) = [abs(x[1]), -0.8*x[2]]

sysG = HD.GeneralSystem(f_ball, h_ball, Δ_ball; direction=-1)
probG1 = HD.prob(sysG, [10.0, 0.0], (0.0, 2.0))
probGz = HD.prob(sysG, [1.0, 0.0], (0.0, 5.0))

@testset "General" begin
    
    solG1 = HD.solve(probG1, HD.RK4())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solG1(1) ≈ [5.095, -9.81] atol=soltol
    # After bounce, position
    @test solG1(1.75)[1] ≈ -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)) atol=soltol
    # One bounce
    @test length(solG1.event_times) == 1
    # At the right time
    @test solG1.event_times[1] ≈ sqrt(10*2/9.81) atol=soltol

    # Throws warning for repeated event
    @test_logs (:warn,) HD.solve(probGz, HD.RK4())

end

