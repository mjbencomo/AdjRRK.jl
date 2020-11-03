# Looking at the accuracy of the forward and adjoint algorithms of RK, IDT, and
# RRK on a skew-symmetric problem that is energy conservative.

# Results:
# - Proper convergence of forward, linearized, and adjoint methods, where true
#   final value is used as the final condition for the adjoint.
# - When using final value computed by respective method as the final condition
#   for the adjoint methods we observed better convergence of RK (one order
#   higher than expected) and time reversibility of the relaxed methods.
# - When using the final value of the linearized methods as the final condition
#   for the adjoint methods we observed slighlty degraded time reversibility
#   results for the relaxed methods.

using AdjRRK
using LinearAlgebra
using Plots
using UnPack

rk = rk4
arrks = AdjRRK_struct()

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
    return 0.5*norm(u)^2
end
function ∇η(u)
    return u
end
function Hη(u,δu;adj=false)
    return δu
end

@pack! arrks = f,df,η,∇η,Hη

ts = Time_struct()
t0  = 0
T   = 100
@pack! ts = t0,T

u0 = [2,1]
u0_lin = u0
uT_adj = exp(A.*T)*u0
@pack! arrks = u0,u0_lin,uT_adj

# Time = (t0,T,dt0/2)
# u_RK,Δη_RK,t_RK = RK_solver(u0,(f,η),Time,rk; return_Δη=true,return_time=true)
# u_IDT,γ_IDT,Δη_IDT = IDT_solver(u0,(f,η,∇η),Time,rk; return_Δη=true)
# u_RRK,γ_RRK,Δη_RRK,t_RRK,dt_corr = RRK_solver(u0,(f,η,∇η),Time,rk; return_Δη=true,return_time=true)
#
# w_RK  = RK_solver((u0,u_RK),(f,df),Time,rk; lin=true)
# w_IDT = IDT_solver((u0,u_IDT,γ_IDT),(f,df,∇η,Hη),Time,rk; lin=true)
# w_RRK = RRK_solver((u0,u_RRK,γ_RRK,dt_corr),(f,df,∇η,Hη),Time,rk; lin=true)
#
# z_RK  = RK_solver((uT,u_RK),(f,df),Time,rk; adj=true)
# z_IDT = IDT_solver((uT,u_IDT,γ_IDT),(f,df,∇η,Hη),Time,rk; adj=true)
# z_RRK = RRK_solver((uT,u_RRK,γ_RRK,dt_corr),(f,df,∇η,Hη),Time,rk; adj=true)
#
# #plotting entropy
# plot(t_RK,Δη_RK,label="RK",marker=:circle)
# plot!(t_RK,Δη_IDT,label="IDT",marker=:square)
# plot!(t_RRK,Δη_RRK,label="RRK",marker=:star)
# display(plot!())
#
# #plotting relaxation parameter
# plot(t_RK,γ_IDT,label="IDT",marker=:square)
# plot!(t_RRK,γ_RRK,label="RRK",marker=:star)
# display(plot!())
#
# #plotting forward solution
# plot(t_RK,u_RK[1,:],label="RK",marker=:circle)
# plot!(t_RK,u_IDT[1,:],label="IDT",marker=:square)
# plot!(t_RRK,u_RRK[1,:],label="RRK",marker=:star)
# display(plot!())
#
# plot(t_RK,u_RK[2,:],label="RK",marker=:circle)
# plot!(t_RK,u_IDT[2,:],label="IDT",marker=:square)
# plot!(t_RRK,u_RRK[2,:],label="RRK",marker=:star)
# display(plot!())
#
# #plotting linearzed solution
# plot(t_RK,w_RK[1,:],label="RK",marker=:circle)
# plot!(t_RK,w_IDT[1,:],label="IDT",marker=:square)
# plot!(t_RRK,w_RRK[1,:],label="RRK",marker=:star)
# display(plot!())
#
# plot(t_RK,w_RK[2,:],label="RK",marker=:circle)
# plot!(t_RK,w_IDT[2,:],label="IDT",marker=:square)
# plot!(t_RRK,w_RRK[2,:],label="RRK",marker=:star)
# display(plot!())
#
# #plotting adjoint solution
# plot(t_RK,z_RK[1,:],label="RK",marker=:circle)
# plot!(t_RK,z_IDT[1,:],label="IDT",marker=:square)
# plot!(t_RRK,z_RRK[1,:],label="RRK",marker=:star)
# display(plot!())
#
# plot(t_RK,z_RK[2,:],label="RK",marker=:circle)
# plot!(t_RK,z_IDT[2,:],label="IDT",marker=:square)
# plot!(t_RRK,z_RRK[2,:],label="RRK",marker=:star)
# display(plot!())

dt0 = 1
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

## RK convergence run
solver! = RK_solver!

err_fwd_RK = zeros(Nref+1)
err_lin_RK = zeros(Nref+1)
err_adj_RK = zeros(Nref+1)

for n=1:Nref+1
    ts.dt = dt[n]

    solver!(arrks,ts,rk)
    solver!(arrks,ts,rk;lin=true)
    solver!(arrks,ts,rk;adj=true)
    @unpack u,u_lin,u_adj = arrks

    err_fwd_RK[n] = norm(u[:,end] - uT_adj)
    err_lin_RK[n] = norm(u_lin[:,end] - uT_adj)
    err_adj_RK[n] = norm(u_adj[:,1] - u0)
end

rates_fwd_RK = comp_rates(err_fwd_RK)
rates_lin_RK = comp_rates(err_lin_RK)
rates_adj_RK = comp_rates(err_adj_RK)


## IDT convergence run
solver = IDT_solver!

err_fwd_IDT   = zeros(Nref+1)
err_lin_IDTγ0 = zeros(Nref+1)
err_lin_IDT   = zeros(Nref+1)
err_adj_IDTγ0 = zeros(Nref+1)
err_adj_IDT   = zeros(Nref+1)

for n=1:Nref+1
    ts.dt = dt[n]

    solver!(arrks,ts,rk)
    @unpack u = arrks
    err_fwd_IDT[n] = norm(u[:,end] - uT_adj)

    arrks.γ_cnst = true
    solver!(arrks,ts,rk;lin=true)
    solver!(arrks,ts,rk;adj=true)
    @unpack u_lin,u_adj = arrks
    err_lin_IDTγ0[n] = norm(u_lin[:,end] - uT_adj)
    err_adj_IDTγ0[n] = norm(u_adj[:,1] - u0)

    arrks.γ_cnst = false
    solver!(arrks,ts,rk;lin=true)
    solver!(arrks,ts,rk;adj=true)
    @unpack u_lin,u_adj = arrks
    err_lin_IDT[n] = norm(u_lin[:,end] - uT_adj)
    err_adj_IDT[n] = norm(u_adj[:,1] - u0)
end

rates_fwd_IDT   = comp_rates(err_fwd_IDT)
rates_lin_IDTγ0 = comp_rates(err_lin_IDTγ0)
rates_lin_IDT   = comp_rates(err_lin_IDT)
rates_adj_IDTγ0 = comp_rates(err_adj_IDTγ0)
rates_adj_IDT   = comp_rates(err_adj_IDT)


## RRK convergence run
solver! = RRK_solver!

err_fwd_RRK   = zeros(Nref+1)
err_lin_RRKγ0 = zeros(Nref+1)
err_lin_RRK   = zeros(Nref+1)
err_adj_RRKγ0 = zeros(Nref+1)
err_adj_RRK   = zeros(Nref+1)

for n=1:Nref+1
    ts.dt = dt[n]

    solver!(arrks,ts,rk)
    @unpack u = arrks
    err_fwd_RRK[n] = norm(u[:,end] - uT_adj)

    arrks.γ_cnst = true
    solver!(arrks,ts,rk;lin=true)
    solver!(arrks,ts,rk;adj=true)
    @unpack u_lin,u_adj = arrks
    err_lin_RRKγ0[n] = norm(u_lin[:,end] - uT_adj)
    err_adj_RRKγ0[n] = norm(u_adj[:,1] - u0)

    arrks.γ_cnst = false
    solver!(arrks,ts,rk;lin=true)
    solver!(arrks,ts,rk;adj=true)
    @unpack u_lin,u_adj = arrks
    err_lin_RRK[n] = norm(u_lin[:,end] - uT_adj)
    err_adj_RRK[n] = norm(u_adj[:,1] - u0)
end

rates_fwd_RRK   = comp_rates(err_fwd_RRK)
rates_lin_RRKγ0 = comp_rates(err_lin_RRKγ0)
rates_lin_RRK   = comp_rates(err_lin_RRK)
rates_adj_RRKγ0 = comp_rates(err_adj_RRKγ0)
rates_adj_RRK   = comp_rates(err_adj_RRK)


## Plots
labels = ["RK" "IDT" "RRK"]
markers = [:circle :square :star]

plot(dt,[err_fwd_RK,err_fwd_IDT,err_fwd_RRK],
    title="Forward convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_fwd_RK,rates_fwd_IDT,rates_fwd_RRK],
    label=labels,
    marker=markers,
    xaxis=:flip,
    title="Forward convergence slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())


labels = ["RK" "IDTγ0" "IDT" "RRKγ0" "RRK"]
markers = [:circle :square :square :star :star]
linestyles = [:solid :solid :dash :solid :dash]

plot(dt,[err_lin_RK,err_lin_IDTγ0,err_lin_IDT,err_lin_RRKγ0,err_lin_RRK],
    title="Linearized convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_lin_RK,rates_lin_IDTγ0,rates_lin_IDT,rates_lin_RRKγ0,rates_lin_RRK],
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xaxis=:flip,
    title="Linearized convergence slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())

plot(dt,[err_adj_RK,err_adj_IDTγ0,err_adj_IDT,err_adj_RRKγ0,err_adj_RRK],
    title="Adjoint convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_adj_RK,rates_adj_IDTγ0,rates_adj_IDT,rates_adj_RRKγ0,rates_adj_RRK],
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xaxis=:flip,
    title="Adjoint convergence slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())

## Time reversal results

# The following results demonstrate the time reversibility of RRK and IDT as well
# as the super convergence behavior of adjoint RK4.

err_RK    = zeros(Nref+1)
err_IDTγ0 = zeros(Nref+1)
err_IDT   = zeros(Nref+1)
err_RRKγ0 = zeros(Nref+1)
err_RRK   = zeros(Nref+1)

for n=1:Nref+1
    ts.dt = dt[n]

    RK_solver!(arrks,ts,rk)
    @unpack u = arrks
    arrks.uT_adj = u[:,end]

    RK_solver!(arrks,ts,rk;adj=true)
    @unpack u_adj = arrks
    err_RK[n] = norm(u_adj[:,1] - u0)

    #IDT
    IDT_solver!(arrks,ts,rk)
    @unpack u = arrks
    arrks.uT_adj = u[:,end]

    arrks.γ_cnst = true
    IDT_solver!(arrks,ts,rk;adj=true)
    @unpack u_adj = arrks
    err_IDTγ0[n] = norm(u_adj[:,1] - u0)

    arrks.γ_cnst = false
    IDT_solver!(arrks,ts,rk;adj=true)
    @unpack u_adj = arrks
    err_IDT[n] = norm(u_adj[:,1] - u0)

    #RRK
    RRK_solver!(arrks,ts,rk)
    @unpack u = arrks
    arrks.uT_adj = u[:,end]

    arrks.γ_cnst = true
    RRK_solver!(arrks,ts,rk;adj=true)
    @unpack u_adj = arrks
    err_RRKγ0[n] = norm(u_adj[:,1] - u0)

    arrks.γ_cnst = false
    RRK_solver!(arrks,ts,rk;adj=true)
    @unpack u_adj = arrks
    err_RRK[n] = norm(u_adj[:,1] - u0)
end

rates_RK    = comp_rates(err_RK)
rates_IDTγ0 = comp_rates(err_IDTγ0)
rates_IDT   = comp_rates(err_IDT)
rates_RRKγ0 = comp_rates(err_RRKγ0)
rates_RRK   = comp_rates(err_RRK)

plot(dt,[err_RK,err_IDTγ0,err_IDT,err_RRKγ0,err_RRK],
    title="Adjoint convergence plot",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xlabel="dt",
    ylabel="Error at t=$T")
display(plot!())

plot([rates_RK,rates_IDTγ0,rates_IDT,rates_RRKγ0,rates_RRK],
    label=labels,
    marker=markers,
    linestyle=linestyles,
    xaxis=:flip,
    title="Adjoint convergence slopes/rates",
    ylabel="Rates/slopes",
    xlabel="refinement index")
display(plot!())
