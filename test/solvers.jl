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

sol1 = [5.095, -9.81]
sol2 = -9.81/2*(1.75-sqrt(10*2/9.81))^2 + 9.81*0.8*sqrt(10*2/9.81)*(1.75-sqrt(10*2/9.81))
solevent = sqrt(10*2/9.81)

@testset "FixedLMMs" begin
    
    solG_ab2 = HD.solve(probG, HD.AdamsBashforth2())

    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solG_ab2(1) ≈ sol1 atol=soltol
    # After bounce, position
    @test solG_ab2(1.75)[1] ≈ sol2 atol=soltol
    # One bounce
    @test length(solG_ab2.event_times) == 1
    # At the right time
    @test solG_ab2.event_times[1] ≈ solevent atol=soltol

    solG_ab3 = HD.solve(probG, HD.AdamsBashforth3())
    @test solG_ab3(1) ≈ sol1 atol=soltol
    @test solG_ab3(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_ab3.event_times) == 1
    @test solG_ab3.event_times[1] ≈ solevent atol=soltol

    solG_BDF2 = HD.solve(probG, HD.BDF2())
    @test solG_BDF2(1) ≈ sol1 atol=soltol
    @test solG_BDF2(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_BDF2.event_times) == 1
    @test solG_BDF2.event_times[1] ≈ solevent atol=soltol

end

@testset "AdaptiveLMMs" begin

    solG_AdaptiveABM2 = HD.solve(probG, HD.AdaptiveABM2())
    @test solG_AdaptiveABM2(1) ≈ sol1 atol=soltol
    @test solG_AdaptiveABM2(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_AdaptiveABM2.event_times) == 1
    @test solG_AdaptiveABM2.event_times[1] ≈ solevent atol=soltol

    solG_AdaptiveABM3 = HD.solve(probG, HD.AdaptiveABM3())  
    @test solG_AdaptiveABM3(1) ≈ sol1 atol=soltol
    @test solG_AdaptiveABM3(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_AdaptiveABM3.event_times) == 1
    @test solG_AdaptiveABM3.event_times[1] ≈ solevent atol=soltol

end

@testset "FixedRK" begin

    solG_BackwardEuler = HD.solve(probG, HD.BackwardEuler())
    @test solG_BackwardEuler(1) ≈ sol1 atol=soltol
    @test solG_BackwardEuler(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_BackwardEuler.event_times) == 1
    @test solG_BackwardEuler.event_times[1] ≈ solevent atol=soltol

    solG_ForwardEuler = HD.solve(probG, HD.ForwardEuler())
    @test solG_ForwardEuler(1) ≈ sol1 atol=soltol
    @test solG_ForwardEuler(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_ForwardEuler.event_times) == 1
    @test solG_ForwardEuler.event_times[1] ≈ solevent atol=soltol

    solG_ModifiedMidpoint = HD.solve(probG, HD.ModifiedMidpoint())
    @test solG_ModifiedMidpoint(1) ≈ sol1 atol=soltol
    @test solG_ModifiedMidpoint(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_ModifiedMidpoint.event_times) == 1
    @test solG_ModifiedMidpoint.event_times[1] ≈ solevent atol=soltol

    solG_ModifiedTrap = HD.solve(probG, HD.ModifiedTrap())
    @test solG_ModifiedTrap(1) ≈ sol1 atol=soltol
    @test solG_ModifiedTrap(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_ModifiedTrap.event_times) == 1
    @test solG_ModifiedTrap.event_times[1] ≈ solevent atol=soltol

    # RK4 is tested in the general.jl test

    solG_RichardsonExtrapolation = HD.solve(probG, HD.RichardsonExtrapolation())
    @test solG_RichardsonExtrapolation(1) ≈ sol1 atol=soltol
    @test solG_RichardsonExtrapolation(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_RichardsonExtrapolation.event_times) == 1
    @test solG_RichardsonExtrapolation.event_times[1] ≈ solevent atol=soltol

    solG_BackwardEuler = HD.solve(probG, HD.BackwardEuler())
    @test solG_BackwardEuler(1) ≈ sol1 atol=soltol
    @test solG_BackwardEuler(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_BackwardEuler.event_times) == 1
    @test solG_BackwardEuler.event_times[1] ≈ solevent atol=soltol

    solG_ImplicitTrap = HD.solve(probG, HD.ImplicitTrap())
    @test solG_ImplicitTrap(1) ≈ sol1 atol=soltol
    @test solG_ImplicitTrap(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_ImplicitTrap.event_times) == 1
    @test solG_ImplicitTrap.event_times[1] ≈ solevent atol=soltol

    solG_RadauIIA = HD.solve(probG, HD.RadauIIA())
    @test solG_RadauIIA(1) ≈ sol1 atol=soltol
    @test solG_RadauIIA(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_RadauIIA.event_times) == 1
    @test solG_RadauIIA.event_times[1] ≈ solevent atol=soltol

end

@testset "AdaptiveRK" begin

    solG_RK23 = HD.solve(probG, HD.RK23())
    @test solG_RK23(1) ≈ sol1 atol=soltol
    @test solG_RK23(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_RK23.event_times) == 1
    @test solG_RK23.event_times[1] ≈ solevent atol=soltol

    solG_RK45 = HD.solve(probG, HD.RK45())
    @test solG_RK45(1) ≈ sol1 atol=soltol
    @test solG_RK45(1.75)[1] ≈ sol2 atol=soltol
    @test length(solG_RK45.event_times) == 1
    @test solG_RK45.event_times[1] ≈ solevent atol=soltol

end