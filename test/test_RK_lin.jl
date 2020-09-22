# Testing linearized and adjoint RK algorithms on a simple system of ODEs
using AdjRRK
using Test
using LinearAlgebra

rk = rk4
C = rand(2,2)
Time = (0,1,0.01)

function f(u)
    return C*u
end

function df(u,δu;adj=false)
    if adj
        return transpose(C)*δu
    end
    return C*δu
end

# Solving ODE system
u0 = randn(2)
u,t = RK_solver(u0,f,Time,rk;return_time=true)

# DERIVATIVE TEST
# Given that the problem is linear, it suffices to compare the forward algorithm
# with its linearization.
@testset "derivative test" begin
    w0 = u[:,1]
    w = RK_solver((w0,u),(f,df),Time,rk;lin=true)
    norm_diff = norm(u-w)/norm(u)
    @test norm_diff<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin
    w0 = randn(2)
    w = RK_solver((w0,u),(f,df),Time,rk;lin=true)
    zT = randn(2)
    z = RK_solver((zT,u),(f,df),Time,rk;adj=true)
    ipt = w[:,1]⋅z[:,1] - w[:,end]⋅z[:,end]
    ipt = abs(ipt)/norm(w)*norm(z[:,end])
    @test ipt<AdjRRK.IPT_TOL
end
