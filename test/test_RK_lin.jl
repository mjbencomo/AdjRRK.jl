# Testing RK algorithm (linearized and adjoint) on a simple system of ODEs.

using AdjRRK
using LinearAlgebra
using Test
using UnPack


# setting RKs and AdjRRK_struct structs
rk = rk4
arrks = AdjRRK_struct()

# setting time stuff
ts = Time_struct()
ts.t0 = 0
ts.T  = 1
ts.dt = 0.01

# setting RHS function and its Jacobian
C = randn(2,2)
function f(u)
    return C*u
end
function df(u,δu;adj=false)
    if adj
        return transpose(C)*δu
    end
    return C*δu
end
@pack! arrks = f,df

# Running RK
arrks.u0 = randn(2)
RK_solver!(arrks,ts,rk)

# DERIVATIVE TEST
# Given that the problem is linear, it suffices to compare the forward algorithm
# with its linearization.
@testset "derivative test" begin
    arrks.u0_lin = arrks.u0
    RK_solver!(arrks,ts,rk;lin=true)
    @unpack u,u_lin = arrks

    norm_diff = norm(u-u_lin)/norm(u)
    @test norm_diff<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin
    arrks.u0_lin = randn(2)
    RK_solver!(arrks,ts,rk;lin=true)
    arrks.uT_adj = randn(2)
    RK_solver!(arrks,ts,rk;adj=true)

    @unpack u_lin,u_adj= arrks
    ipt = u_lin[:,1]⋅u_adj[:,1] - u_lin[:,end]⋅u_adj[:,end]
    ipt = abs(ipt)
    ipt /= norm(u_lin)*norm(u_adj[:,end])
    @test ipt<AdjRRK.IPT_TOL
end
