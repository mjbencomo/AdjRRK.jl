# Testing convergence of RK, RRK, and IDT methods on conserved exponential
# entropy example, see Ranocha et al 2019.

using AdjRRK
using LinearAlgebra
using Test

function avg(x)
    return sum(x)/length(x)
end

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

function true_sol(t)
    u = zeros(2)
    e = exp(1)
    E = √e + e

    u[2] = log( (exp(E*t)*E)/(√e + exp(E*t)) )
    u[1] = log( e+exp(3/2) ) - log( √e + exp(E*t) )

    return u
end

t0 = 0
T  = 5
dt0 = 0.5
u0 = [1,0.5]
u_true = true_sol(T)

Nref = 4
dt = zeros(Nref+1)
dt[1] = dt0
for n=1:Nref
    dt[n+1] = dt[n]/2
end


## RK convergence run
solver = RK_solver
ops = f

#RK4
rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true)
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_RK4 = err
rates_RK4 = rates

#RK2 = Heun's method
rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true)
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_RK2 = err
rates_RK2 = rates

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
solver = IDT_solver
ops = (f,η,∇η)

#RK4
rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true[:])
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_IDT4 = err
rates_IDT4 = rates

#RK2 = Heun's method
rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true[:])
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_IDT2 = err
rates_IDT2 = rates

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
solver = RRK_solver
ops = (f,η,∇η)

#RK4
rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true[:])
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_RRK4 = err
rates_RRK4 = rates

#RK2 = Heun's method
rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)
    err[n] = norm(u[:,end] - u_true[:])
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_RRK2 = err
rates_RRK2 = rates

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
