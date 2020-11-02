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

# Running RK
u0 = [1.5,1]
@pack! arrks = u0
RK_solver!(arrks,ts,rk)
@unpack u = arrks

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
    arrks.u0_lin = d
    RK_solver!(arrks,ts,rk;lin=true)
    @unpack u_lin = arrks

    arrks_h = AdjRRK_struct()
    arrks_h.f  = arrks.f
    arrks_h.df = arrks.df
    for n=1:Nref+1
        arrks_h.u0 = u0 + h[n].*d
        RK_solver!(arrks_h,ts,rk)
        uh = arrks_h.u

        FD_err[n] = norm((uh-u)./h[n] - u_lin)
        if n<Nref+1
            h[n+1] = h[n]/2
        end


        # uh0 = u[:,1] + h[n].*d
        # uh  = RK_solver!(uh0,f,Time,rk)
        # FD_err[n] = norm((uh-u)./h[n] - w)
        # if n<Nref+1
        #     h[n+1] = h[n]/2
        # end
    end
    rate = log2(FD_err[Nref])-log2(FD_err[Nref+1])
    @test abs(rate-1)<AdjRRK.DRV_TOL
end

# INNER PRODUCT TEST
@testset "inner product test" begin
    # Running linear RK
    arrks.u0_lin = randn(2)
    RK_solver!(arrks,ts,rk;lin=true)

    # Running adjoint RK
    arrks.uT_adj = randn(2)
    RK_solver!(arrks,ts,rk;adj=true)

    @unpack u_lin,u_adj = arrks
    ipt = u_lin[:,1]⋅u_adj[:,1] - u_lin[:,end]⋅u_adj[:,end]
    ipt = abs(ipt)
    ipt /= norm(u_lin)*norm(u_adj[:,end])
    @test ipt<AdjRRK.IPT_TOL
end
