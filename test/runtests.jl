using SafeTestsets
@safetestset "Testing adjoint RK on linear ODE" begin
        include("test_RK_lin.jl")
end
@safetestset "Testing adjoint RK on nonlinear pendulum" begin
        include("test_RK_nlinpen.jl")
end
@safetestset "Testing adjoint IDT on nonlinear pendulum" begin
        include("test_IDT_nlinpen.jl")
end
@safetestset "Testing adjoint RRK on nonlinear pendulum" begin
        include("test_RRK_nlinpen.jl")
end
@safetestset "Testing convergence on conserved exp entropy problem" begin
        include("test_RRK_nlinpen.jl")
end
