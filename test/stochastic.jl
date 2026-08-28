# # Comment out when not testing within the file to avoid redundant imports
import HybridDynamics as HD
using Test

# Define problem with known solution
function f_ball(x, t)
    g = 9.81
    q, v = x
    return [v, -g]
end
w(x, t) = [0.0; 0.02;;]
h_ball(x) = x[1]
Δ_ball(x) = [abs(x[1]), -0.8*x[2]]

sysS = HD.StochasticSystem(f_ball, w, h_ball, Δ_ball; direction=-1)
probS = HD.prob(sysS, [10.0, 0.0], (0.0, 2.0))

@testset "Stochastic" begin
    
    solS = HD.solve(probS)

    # Test to verify the solution at t=1 is around an exact value
    @test solS.x[101] ≈ [5.095, -9.81] atol=0.5
    # After bounce, state
    @test solS.x[177] ≈ [-9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81)), 8.1] atol=2
    # One bounce
    @test length(solS.event_times) == 1
    # At the right time
    @test solS.event_times[1] ≈ sqrt(10*2/9.81) atol=1
    # Don't allow interpolation
    @test_throws ErrorException solS(1.0)

end

