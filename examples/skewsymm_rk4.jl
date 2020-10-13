# Testing the accuracy of RK4, and it's linearization and adjoint, on a simple
# linear skew-symmetric problem.

# Convergence of order 5 is observed for the adjoint problem whenever the computed
# solution for the forward problem is used as the final time condition.
# However, convergence of order 4 is observed if the true final condition is used.

using AdjRRK
using LinearAlgebra
using Plots

rk = rk4
A = [0 1; -1 0]


function f(u)
    return A*u
end
function df(u,δu;adj=false)
    if adj
        return transpose(A)*δu
    end
    return A*δu
end

function η(u)
    return 0.5*norm(u)
end
function ∇η(u)
    return u
end
function Hη(u,δu;adj=false)
    return δu
end

function true_sol(T,u0)
    return exp(A.*T)*u0
end

t0  = 0
T   = 100
dt0 = 1

u0 = [2,1]
uT = true_sol(T,u0)


Nref = 6
dt = zeros(Nref+1)
dt[1] = dt0
for n=1:Nref
    dt[n+1] = dt[n]/2
end

function comp_rates(err)
    rates = zeros(Nref)
    for n=1:Nref
        rates[n] = log2(err[n]) - log2(err[n+1])
    end
    return rates
end


err_fwd_RK = zeros(Nref+1)
err_lin_RK = zeros(Nref+1)
err_adj_RK = zeros(Nref+1)

for n=1:Nref+1
    Time = (t0,T,dt[n])
    u = RK_solver(u0,f,Time,rk)
    w = RK_solver((u0,u),(f,df),Time,rk;lin=true)
    z = RK_solver((uT,u),(f,df),Time,rk;adj=true) #yields order-4 convergence
    #z = RK_solver((u[:,end],u),(f,df),Time,rk;adj=true) #yields order-5 convergence

    err_fwd_RK[n] = norm(u[:,end] - uT)
    err_lin_RK[n] = norm(w[:,end] - uT)
    err_adj_RK[n] = norm(z[:,1] - u0)
end

rates_fwd_RK = comp_rates(err_fwd_RK)
rates_lin_RK = comp_rates(err_lin_RK)
rates_adj_RK = comp_rates(err_adj_RK)


# Plots
plot(dt,err_fwd_RK,
    title="RK4: convergence plot",
    xaxis=:log,
    yaxis=:log,
    marker=:circle,
    legend=false,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot(rates_fwd_RK,
    title="RK4: convergence slopes/rates",
    marker=:circle,
    xaxis=:flip,
    legend=false,
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())

plot(dt,err_lin_RK,
    title="Linearized RK4: convergence plot",
    xaxis=:log,
    yaxis=:log,
    marker=:circle,
    legend=false,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot(rates_lin_RK,
    title="Linearized RK4: convergence slopes/rates",
    marker=:circle,
    xaxis=:flip,
    legend=false,
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())

plot(dt,err_adj_RK,
    title="Adjoint RK4: convergence plot",
    xaxis=:log,
    yaxis=:log,
    marker=:circle,
    legend=false,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot(rates_adj_RK,
    title="Adjoint RK4: convergence slopes/rates",
    marker=:circle,
    xaxis=:flip,
    legend=false,
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())
