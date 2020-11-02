# Dissipative exponential entropy example, see Ranocha et al 2019

using AdjRRK
using LinearAlgebra
using Plots, Printf
using UnPack

arrks = AdjRRK_struct()

function f(u)
    return @. -exp(u)
end

function df(u,δu;adj=false)
    J = @. -exp(u)
    return J*δu
end

function η(u)
    return exp(u[1])
end

function ∇η(u)
    return exp.(u)
end

function Hη(u,δu;adj=false)
    H = exp.(u)
    return H*δu
end
@pack! arrks = f,df,η,∇η,Hη


function true_sol(t)
    u = zeros(1)
    u[1] = -log(exp(-1/2)+t)
    return u
end

ts = Time_struct()
t0 = 0
T  = 5
@pack! ts = t0,T

arrks.u0 = [0.5]
u_true = true_sol(T)

Nref = 10
dt0 = 0.5


## RK convergence run
errs_RK4,rates_RK4,dt = AdjRRK.conv_test!(RK_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_RK2,rates_RK2,dt = AdjRRK.conv_test!(RK_solver!,arrks,ts,rk2,dt0,Nref,u_true)

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[errs_RK4,errs_RK2],
    title="RK convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_RK4,rates_RK2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="RK slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())


## IDT convergence run
errs_IDT4,rates_IDT4,dt = AdjRRK.conv_test!(IDT_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_IDT2,rates_IDT2,dt = AdjRRK.conv_test!(IDT_solver!,arrks,ts,rk2,dt0,Nref,u_true)

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[errs_IDT4,errs_IDT2],
    title="IDT convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_IDT4,rates_IDT2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="IDT slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())


## RRK convergence run
errs_RRK4,rates_RRK4,dt = AdjRRK.conv_test!(RRK_solver!,arrks,ts,rk4,dt0,Nref,u_true)
errs_RRK2,rates_RRK2,dt = AdjRRK.conv_test!(RRK_solver!,arrks,ts,rk2,dt0,Nref,u_true)

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[errs_RRK4,errs_RRK2],
    title="RRK convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_RRK4,rates_RRK2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="RRK slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())
