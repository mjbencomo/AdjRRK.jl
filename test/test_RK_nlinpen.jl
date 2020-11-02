# Testing out RK algorithm (linearized and adjoint code) on a small nonlinear ODE
# system (nolinear pendulum) from Ranocha 2019.

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

# initial condition
arrks.u0 = [1.5,1]

@testset "derivative test" begin
    arrks.u0_lin = randn(2)
    Nref = 10
    h0 = 2^(-8)

    arrks_h = AdjRRK_struct(arrks)
    @pack! arrks_h = f,df
    errs,rate,h = AdjRRK.derv_test!(arrks,arrks_h,ts,rk,h0,Nref)
    @test abs(rate[end]-1)<AdjRRK.DRV_TOL
end

@testset "inner product test" begin
    ipt = AdjRRK.ip_test!(arrks,ts,rk)
    @test ipt<AdjRRK.IPT_TOL
end
