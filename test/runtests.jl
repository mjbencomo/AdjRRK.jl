using SafeTestsets

IPT_TOL = 1e-13
DRV_TOL = 1e-1

@safetestset "Testing adjoint RK on linear ODE." begin include("test_RK_lin.jl") end


# using AdjRRK
# using Test
#
#
#
# @testset "AdjRRK.jl" begin
#     @test IPT_TOL<1
# end
