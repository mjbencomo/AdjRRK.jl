# Testing convergence of RK, RRK, and IDT methods on conserved exponential
# entropy example, see Ranocha et al 2019.

using AdjRRK
using LinearAlgebra
using Test
using UnPack

function avg(x)
    return sum(x)/length(x)
end

# initializing AdjRRK_struct
arrks = AdjRRK_struct()

function f(u)
    return [-exp(u[2]),exp(u[1])]
end
function df(u,δu;adj=false)
    J = [0 -exp(u[2]); exp(u[1]) 0]
    if adj
        J = transpose(J)
    end
    return J*δu
end
@pack! arrks = f,df

function η(u)
    return exp(u[1])+exp(u[2])
end
function ∇η(u)
    return [exp(u[1]),exp(u[2])]
end
function Hη(u,δu;adj=false)
    H = [exp(u[1]) 0; 0 exp(u[2])]
    if adj
        H = transpose(H)
    end
    return H*δu
end
@pack! arrks = η,∇η,Hη

function true_sol(t)
    u = zeros(2)
    e = exp(1)
    E = √e + e

    u[2] = log( (exp(E*t)*E)/(√e + exp(E*t)) )
    u[1] = log( e+exp(3/2) ) - log( √e + exp(E*t) )

    return u
end


# setting time stuff
ts = Time_struct()
t0 = 0
T  = 5
@pack! ts = t0,T

dt0 = 0.5
Nref = 4

arrks.u0 = [1,0.5]
u_true = true_sol(T)

## RK convergence run
errs_rk4,rates_RK4,dt = AdjRRK.conv_test!(RK_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_rk2,rates_RK2,dt = AdjRRK.conv_test!(RK_solver!,arrks,ts,rk2,dt0,Nref,u_true)

@testset "RK convergence test" begin
    @test abs(avg(rates_RK4)-4) < AdjRRK.CNV_TOL
    @test abs(avg(rates_RK2)-2) < AdjRRK.CNV_TOL
end

# labels = ["RK4" "Heun"]
# markers = [:circle :square]
# plot(dt,[err_RK4,err_RK2],
#     title="RK convergence plot",
#     xaxis=:log,
#     yaxis=:log,
#     label=labels,
#     marker=markers,
#     xlabel="dt",
#     ylabel="Error at t=$T")
# display(plot!())
#
# plot([rates_RK4,rates_RK2],
#     label=labels,
#     marker=markers,
#     xaxis=:flip,
#     title="RK slopes/rates",
#     ylabel="Rates/slopes",
#     xlabel="refinement index")
# display(plot!())


## IDT convergence run
errs_IDT4,rates_IDT4,dt = AdjRRK.conv_test!(IDT_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_IDT2,rates_IDT2,dt = AdjRRK.conv_test!(IDT_solver!,arrks,ts,rk2,dt0,Nref,u_true)

@testset "IDT convergence test" begin
    @test abs(avg(rates_IDT4)-3) < AdjRRK.CNV_TOL
    @test abs(avg(rates_IDT2)-1) < AdjRRK.CNV_TOL
end

# labels = ["RK4" "Heun"]
# markers = [:circle :square]
# plot(dt,[err_IDT4,err_IDT2],
#     title="IDT convergence plot",
#     xaxis=:log,
#     yaxis=:log,
#     label=labels,
#     marker=markers,
#     xlabel="dt",
#     ylabel="Error at t=$T")
# display(plot!())
#
# plot([rates_IDT4,rates_IDT2],
#     label=labels,
#     marker=markers,
#     xaxis=:flip,
#     title="IDT slopes/rates",
#     ylabel="Rates/slopes",
#     xlabel="refinement index")
# display(plot!())

## RRK convergence run
errs_RRK4,rates_RRK4,dt = AdjRRK.conv_test!(RRK_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_RRK2,rates_RRK2,dt = AdjRRK.conv_test!(RRK_solver!,arrks,ts,rk2,dt0,Nref,u_true)

@testset "RRK convergence test" begin
    @test abs(avg(rates_RRK4)-4) < AdjRRK.CNV_TOL
    @test abs(avg(rates_RRK2)-2) < AdjRRK.CNV_TOL
end

# labels = ["RK4" "Heun"]
# markers = [:circle :square]
# plot(dt,[err_RRK4,err_RRK2],
#     title="RRK convergence plot",
#     xaxis=:log,
#     yaxis=:log,
#     label=labels,
#     marker=markers,
#     xlabel="dt",
#     ylabel="Error at t=$T")
# display(plot!())
#
# plot([rates_RRK4,rates_RRK2],
#     label=labels,
#     marker=markers,
#     xaxis=:flip,
#     title="RRK slopes/rates",
#     ylabel="Rates/slopes",
#     xlabel="refinement index")
# display(plot!())
