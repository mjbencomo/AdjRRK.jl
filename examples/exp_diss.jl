# Dissipative exponential entropy example, see Ranocha et al 2019

using AdjRRK
using Plots, Printf

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

function true_sol(t)
    u = zeros(1)
    u[1] = -log(exp(-1/2)+t)
    return u
end

t0 = 0
T  = 5
dt0 = 0.5
u0 = [0.5]
u_true = true_sol(T)

Nref = 10
dt = zeros(Nref+1)
dt[1] = dt0


## RK convergence run
solver = RK_solver
ops = f

rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true)
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk4 = err
rates_rk4 = rates

rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true)
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk2 = err
rates_rk2 = rates

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[err_rk4,err_rk2],
    title="RK convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_rk4,rates_rk2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="RK slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())


## IDT convergence run

solver = IDT_solver
ops = (f,η,∇η)

rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true[:])
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk4 = err
rates_rk4 = rates

rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true[:])
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk2 = err
rates_rk2 = rates

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[err_rk4,err_rk2],
    title="IDT convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_rk4,rates_rk2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="IDT slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())

## RRK convergence run

solver = RRK_solver
ops = (f,η,∇η)

rk = rk4
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true[:])
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk4 = err
rates_rk4 = rates

rk = rk2
err = zeros(Nref+1)
for n=1:Nref+1
    Time = (t0,T,dt[n])
    u,γ = solver(u0,ops,Time,rk)

    err[n] = norm(u[:,end] - u_true[:])
    if n<Nref+1
        dt[n+1] = dt[n]/2
    end
end
rates = zeros(Nref)
for n=1:Nref
    rates[n] = log2(err[n]) - log2(err[n+1])
end
err_rk2 = err
rates_rk2 = rates

labels = ["RK4" "Heun"]
markers = [:circle :square]
plot(dt,[err_rk4,err_rk2],
    title="RRK convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_rk4,rates_rk2],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="RRK slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())
