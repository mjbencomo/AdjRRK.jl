# Testing out IDT algorithm (linearized and adjoint code) on a small nonlinear ODE
# system (nlinear pendulum) from Ranocha 2019.

using AdjRRK
using LinearAlgebra
using Test
using UnPack

# setting RKs and AdjRRK_struct structs
rk = rk4
arrks = AdjRRK_struct()

# setting time stuff
ts = Time_struct()
t0 = 0
T  = 500
dt = 0.9
@pack! ts = t0,T,dt

# setting RHS function and its Jacobian
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
@pack! arrks = f,df

# setting entropy functions
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
@pack! arrks = η,∇η,Hη

#initial condition
arrks.u0 = [1.5,1]

@testset "entropy stability test" begin
    arrks.return_Δη=true
    IDT_solver!(arrks,ts,rk)
    @unpack Δη = arrks
    @test abs(Δη[end])<AdjRRK.ENT_TOL
end

@testset "derivative test" begin
    arrks.u0_lin = randn(2)
    Nref = 8
    h0 = 2^(-6)

    # # γcnst=true case
    # arrks.γcnst = true
    # arrks_h = AdjRRK_struct(arrks)
    # @pack! arrks_h = f,df,η,∇η,Hη
    # errs_γ0,rate_γ0,h = Addj.derv_test!(arrks,arrks_h,ts,rk,h0,Nref)
    # @test abs(rate_γ0[end]-1)<AdjRRK.DRV_TOL

    # γcnst=false case
    arrks.γcnst = false
    arrks_h = AdjRRK_struct(arrks)
    @pack! arrks_h = f,df,η,∇η,Hη
    errs,rate,h = AdjRRK.derv_test!(arrks,arrks_h,ts,rk,h0,Nref)
    @test abs(rate[end]-1)<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin

    # γcnst=true case
    arrks.γcnst = true
    ipt_γ0 = AdjRRK.ip_test!(arrks,ts,rk)
    @test ipt_γ0<AdjRRK.IPT_TOL

    # γcnst=false case
    arrks.γcnst = false
    ipt = AdjRRK.ip_test!(arrks,ts,rk)
    @test ipt<AdjRRK.IPT_TOL
end
