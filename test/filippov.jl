# # Comment out when not testing within the file to avoid redundant imports
# import HybridDynamics as HD
# using Test

# # Tolerance for testing the output of solution structs
# soltol = 1e-2

@testset "Filippov" begin

    # Define problem with known solution
    F(x) = [1, x[1]]
    G(x) = [0, 1]
    H(x) = x[2]

    sysF = HD.FilippovSystem(F, G, H)
    probF = HD.prob(sysF, [-1.5, 1.0], (0.0, 3.0))
    solF = HD.solve(probF, HD.RK4())
    
    # Test to verify the solution at t=1 is within a tolerance of an exact value
    @test solF(1) ≈ [-0.5, 0] atol=soltol
    # After slide, position
    @test solF(2.625) ≈ [1, 0.5] atol=soltol

    # Throws sliding info
    @test_logs (:info, r"Sliding mode entered at t = .*") (:info, r"Sliding mode exited at t = .*") HD.solve(probF, HD.RK4(); track_sliding=:both)

end
