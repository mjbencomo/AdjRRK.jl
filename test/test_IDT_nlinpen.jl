# Testing out IDT algorithm (linearized and adjoint code) on a small nonlinear ODE
# system (nlinear pendulum) from Ranocha 2019.

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

function η(u)
    return 0.5*u[1]^2 - cos(u[2])
end

function ∇η(u)
    d1 = u[1]
    d2 = sin(u[2])
    return [d1,d2]
end

function Hη(u,δu;adj=false)
    H = [1 0; 0 cos(u[2])]
    if adj
        H = transpose(H)
    end
    return H*δu
end

# Running RRK
u0 = [1.5,1]
u,γ,Δηu,t = IDT_solver(u0,(f,η,∇η),Time,rk;
            return_time=true,
            return_Δη=true)

# ENTROPY STABILITY TEST
@testset "entropy stability test" begin
    @test abs(Δηu[end])<AdjRRK.ENT_TOL
end

# DERIVATIVE TEST
@testset "derivative test" begin
    d = randn(2)
    Nref = 15

    # This test is commented out since the derivative test is expected to fail.
    # # IDTγ0
    # h = 2^(-6)
    # FD_err = zeros(Nref+1)
    # w = IDT_solver((d,u,γ),(f,df),Time,rk;lin=true,γcnst=true)
    # for n=1:Nref+1
    #     uh0 = u[:,1] + h.*d
    #     uh,γh = IDT_solver(uh0,(f,η,∇η),Time,rk;γcnst=true)
    #     FD_err[n] = norm((uh[:,end]-u[:,end])./h - w[:,end])
    #     h = h/2
    # end
    # rate_γ0 = log2(FD_err[Nref]) - log2(FD_err[Nref+1])
    # @test abs(rate_γ0-1)<AdjRRK.DRV_TOL

    # IDT
    h = 2^(-6)
    FD_err = zeros(Nref+1)
    w = IDT_solver((d,u,γ),(f,df,∇η,Hη),Time,rk;lin=true)
    for n=1:Nref+1
        uh0 = u[:,1] + h.*d
        uh,γh = IDT_solver(uh0,(f,η,∇η),Time,rk)
        FD_err[n] = norm((uh[:,end]-u[:,end])./h - w[:,end])
        h = h/2
    end
    rate = log2(FD_err[Nref]) - log2(FD_err[Nref+1])
    @test abs(rate-1)<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin

    w0 = randn(2)
    zT = randn(2)

    # IDTγ0
    w = IDT_solver((w0,u,γ),(f,df),Time,rk;lin=true,γcnst=true)
    z = IDT_solver((zT,u,γ),(f,df),Time,rk;adj=true,γcnst=true)
    ipt_γ0 = w[:,1]⋅z[:,1] - w[:,end]⋅z[:,end]
    ipt_γ0 = abs(ipt_γ0)/norm(w)*norm(z[:,end])
    @test ipt_γ0 < AdjRRK.IPT_TOL

    # IDT
    w = IDT_solver((w0,u,γ),(f,df,∇η,Hη),Time,rk;lin=true)
    z = IDT_solver((zT,u,γ),(f,df,∇η,Hη),Time,rk;adj=true)
    ipt = w[:,1]⋅z[:,1] - w[:,end]⋅z[:,end]
    ipt = abs(ipt)/norm(w)*norm(z[:,end])
    @test ipt < AdjRRK.IPT_TOL
end
