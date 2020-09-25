using SafeTestsets
@safetestset "Testing adjoint RK on linear ODE" begin
        include("test_RK_lin.jl")
end
@safetestset "Testing adjoint RK on nonlinear pendulum ODE" begin
        include("test_RK_nlinpen.jl")
end
@safetestset "Testing adjoint IDT on nonlinear pendulum ODE" begin
        include("test_IDT_nlinpen.jl")
end
