# Testing out RK algorithm (linearized and adjoint code) on a small nonlinear ODE
# system (nolinear pendulum) from Ranocha 2019.

using AdjRRK
using LinearAlgebra
using Test

rk = rk4
Time = (0,500,0.9)

function f(u)
    return [-sin(u[2]),u[1]]
end

function df(u,δu;adj=false)
    J = [0 -cos(u[2]); 1 0]
    if adj
        J = transpose(J)
    end
    return J*δu
end

# Running RK
u0 = [1.5,1]
u = RK_solver(u0,f,Time,rk)

# DERIVATIVE TEST
# Given that the problem is nonlinear, we compare linearization output with a
# finite difference approximation and compute the convergence rate, which should
# be one.
@testset "derivative test" begin
    d      = randn(2)
    Nref   = 10
    FD_err = zeros(Nref+1)
    h      = zeros(Nref+1)
    h[1]   = 2^(-8)
    w = RK_solver((d,u),(f,df),Time,rk;lin=true)

    for n=1:Nref+1
        uh0 = u[:,1] + h[n].*d
        uh  = RK_solver(uh0,f,Time,rk)
        FD_err[n] = norm((uh-u)./h[n] - w)
        if n<Nref+1
            h[n+1] = h[n]/2
        end
    end
    rate = log2(FD_err[Nref])-log2(FD_err[Nref+1])
    @test abs(rate-1)<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin
    # Running linear RK
    w0 = randn(2)
    w  = RK_solver((w0,u),(f,df),Time,rk;lin=true)

    # Running adjoint RK
    zT = randn(2)
    z = RK_solver((zT,u),(f,df),Time,rk;adj=true)

    ipt = w[:,1]'*z[:,1] - w[:,end]'*z[:,end]
    ipt = abs(ipt)/norm(w)*norm(z[:,end])
    @test ipt<AdjRRK.IPT_TOL
end
