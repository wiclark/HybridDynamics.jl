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
probG = HD.prob(sysG, [10.0, 0.0], (0.0, 2.0))

@testset "General" begin
    
    solG = HD.solve(probG, HD.RK4())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solG(1) ≈ [5.095, -9.81] atol=soltol
    # After bounce, position
    @test solG(1.75)[1] ≈ -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)) atol=soltol
    # One bounce
    @test length(solG.event_times) == 1
    # At the right time
    @test solG.event_times[1] ≈ sqrt(10*2/9.81) atol=soltol

end

