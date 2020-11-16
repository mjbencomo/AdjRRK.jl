# Testing out RRK method (linearized and adjoint code) on a small nonlinear ODE
# system (nlinear pendulum) from Ranocha 2019.

# Comparing the linearized and adjoint solutions from RRK algs when relaxation
# parameter γ is viewed as a constant or function of state fields during
# linearization. Also, looking at the effects of considering the last time step
# correction in the linearization.

using AdjRRK
using LinearAlgebra
using Plots
using UnPack

λ=0
function f(u)
    f1 = -sin(u[2]) -λ*u[1]
    f2 = u[1]
    return [f1,f2]
end
function df(u,δu;adj=false)
    J = [-λ -cos(u[2]); 1 0]
    if adj
        return transpose(J)*δu
    end
    return J*δu
end

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
        return transpose(H)*δu
    end
    return H*δu
end

t0 = 0
T  = 500
dt = 0.5
u0 = [1.5,1]

##RK example

arrks_RK = AdjRRK_struct()
@pack! arrks_RK = f,df,η,∇η,Hη
@pack! arrks_RK = u0
arrks_RK.return_time = true
arrks_RK.return_Δη = true

ts_RK = Time_struct()
@pack! ts_RK = t0,T,dt


#running RK
RK_solver!(arrks_RK,ts_RK,rk4)
u1_RK = arrks_RK.u[1,:]
u2_RK = arrks_RK.u[2,:]

plot(u1_RK,u2_RK,
    title="u1 vs u2 (RK)",
    xlabel="u1",
    ylabel="u2",
    legend=false)
display(plot!())

plot(ts_RK.t,arrks_RK.Δη,
    title="Entropy conservation (RK)",
    xlabel="t",
    ylabel="Δη",
    legend=false)
display(plot!())


#running lin RK
arrks_RK.u0_lin = u0
RK_solver!(arrks_RK,ts_RK,rk4;lin=true)
w1_RK = arrks_RK.u_lin[1,:]
w2_RK = arrks_RK.u_lin[2,:]

plot(w1_RK,w2_RK,
    title="w1 vs w2 (lin. RK)",
    xlabel="w1",
    ylabel="w2",
    legend=false)
display(plot!())


#running adj RK
arrks_RK.uT_adj = arrks_RK.u[:,end]
RK_solver!(arrks_RK,ts_RK,rk4;adj=true)
z1_RK = arrks_RK.u_adj[1,:]
z2_RK = arrks_RK.u_adj[2,:]

plot(z1_RK,z2_RK,
    title="z1 vs z2 (adj. RK)",
    xlabel="z1",
    ylabel="z2",
    legend=false)
display(plot!())

##IDT example

arrks_IDT = AdjRRK_struct()
@pack! arrks_IDT = f,df,η,∇η,Hη
@pack! arrks_IDT = u0
arrks_IDT.return_time = true
arrks_IDT.return_Δη = true

ts_IDT = Time_struct()
@pack! ts_IDT = t0,T,dt


#running IDT
IDT_solver!(arrks_IDT,ts_IDT,rk4)
u1_IDT = arrks_IDT.u[1,:]
u2_IDT = arrks_IDT.u[2,:]

plot(u1_IDT,u2_IDT,
    title="u1 vs u2 (IDT)",
    xlabel="u1",
    ylabel="u2",
    legend=false)
display(plot!())

plot(ts_IDT.t,arrks_IDT.Δη,
    title="Entropy conservation (IDT)",
    xlabel="t",
    ylabel="Δη",
    legend=false)
display(plot!())


#running lin IDT (γ0)
arrks_IDT.γ_cnst = true
arrks_IDT.u0_lin = u0
IDT_solver!(arrks_IDT,ts_IDT,rk4;lin=true)
w1_IDTγ0 = arrks_IDT.u_lin[1,:]
w2_IDTγ0 = arrks_IDT.u_lin[2,:]

plot(w1_IDTγ0,w2_IDTγ0,
    title="w1 vs w2 (lin. IDT, γ0)",
    xlabel="w1",
    ylabel="w2",
    legend=false)
display(plot!())


#running lin IDT
arrks_IDT.γ_cnst = false
arrks_IDT.u0_lin = u0
IDT_solver!(arrks_IDT,ts_IDT,rk4;lin=true)
w1_IDT = arrks_IDT.u_lin[1,:]
w2_IDT = arrks_IDT.u_lin[2,:]

plot(w1_IDT,w2_IDT,
    title="w1 vs w2 (lin. IDT)",
    xlabel="w1",
    ylabel="w2",
    legend=false)
display(plot!())


#running adj IDT (γ0)
arrks_IDT.γ_cnst = true
arrks_IDT.uT_adj = arrks_IDT.u[:,end]
IDT_solver!(arrks_IDT,ts_IDT,rk4;adj=true)
z1_IDTγ0 = arrks_IDT.u_adj[1,:]
z2_IDTγ0 = arrks_IDT.u_adj[2,:]

plot(z1_IDTγ0,z2_IDTγ0,
    title="z1 vs z2 (adj. IDT, γ0)",
    xlabel="z1",
    ylabel="z2",
    legend=false)
display(plot!())


#running adj IDT
arrks_IDT.γ_cnst = false
arrks_IDT.uT_adj = arrks_IDT.u[:,end]
IDT_solver!(arrks_IDT,ts_IDT,rk4;adj=true)
z1_IDT = arrks_IDT.u_adj[1,:]
z2_IDT = arrks_IDT.u_adj[2,:]

plot(z1_IDT,z2_IDT,
    title="z1 vs z2 (adj. IDT)",
    xlabel="z1",
    ylabel="z2",
    legend=false)
display(plot!())


##RRK example

arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη
@pack! arrks_RRK = u0
arrks_RRK.return_time = true
arrks_RRK.return_Δη = true

ts_RRK = Time_struct()
@pack! ts_RRK = t0,T,dt


#running RRK
RRK_solver!(arrks_RRK,ts_RRK,rk4)
u1_RRK = arrks_RRK.u[1,:]
u2_RRK = arrks_RRK.u[2,:]

plot(u1_RRK,u2_RRK,
    title="u1 vs u2 (RRK)",
    xlabel="u1",
    ylabel="u2",
    legend=false)
display(plot!())

plot(ts_RRK.t,arrks_RRK.Δη,
    title="Entropy conservation (RRK)",
    xlabel="t",
    ylabel="Δη",
    legend=false)
display(plot!())


#running lin RRK (γ0)
arrks_RRK.γ_cnst = true
arrks_RRK.u0_lin = u0
RRK_solver!(arrks_RRK,ts_RRK,rk4;lin=true)
w1_RRKγ0 = arrks_RRK.u_lin[1,:]
w2_RRKγ0 = arrks_RRK.u_lin[2,:]

plot(w1_RRKγ0,w2_RRKγ0,
    title="w1 vs w2 (lin. RRK, γ0)",
    xlabel="w1",
    ylabel="w2",
    legend=false)
display(plot!())


#running lin RRK
arrks_RRK.γ_cnst = false
arrks_RRK.u0_lin = u0
RRK_solver!(arrks_RRK,ts_RRK,rk4;lin=true)
w1_RRK = arrks_RRK.u_lin[1,:]
w2_RRK = arrks_RRK.u_lin[2,:]

plot(w1_RRK,w2_RRK,
    title="w1 vs w2 (lin. RRK)",
    xlabel="w1",
    ylabel="w2",
    legend=false)
display(plot!())


#running adj RRK (γ0)
arrks_RRK.γ_cnst = true
arrks_RRK.uT_adj = arrks_RRK.u[:,end]
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
z1_RRKγ0 = arrks_RRK.u_adj[1,:]
z2_RRKγ0 = arrks_RRK.u_adj[2,:]

plot(z1_RRKγ0,z2_RRKγ0,
    title="z1 vs z2 (adj. RRK, γ0)",
    xlabel="z1",
    ylabel="z2",
    legend=false)
display(plot!())


#running adj RRK
arrks_RRK.γ_cnst = false
arrks_RRK.uT_adj = arrks_RRK.u[:,end]
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
z1_RRK = arrks_RRK.u_adj[1,:]
z2_RRK = arrks_RRK.u_adj[2,:]

plot(z1_RRK,z2_RRK,
    title="z1 vs z2 (adj. RRK)",
    xlabel="z1",
    ylabel="z2",
    legend=false)
display(plot!())


##derivative tests

Nref = 5
h0 = 2^(-8)

arrks_h = AdjRRK_struct()
@pack! arrks_h = f,df,η,∇η,Hη


##IDT
arrks_IDT.u0_lin = randn(2)

# assuming γ is cnst in linearization
arrks_IDT.γ_cnst = true
arrks_h.γ_cnst = true
errs_IDTγ0,rate_IDTγ0,h = AdjRRK.derv_test!(IDT_solver!,arrks_IDT,arrks_h,ts_IDT,rk4,h0,Nref)

# with proper linearization
arrks_IDT.γ_cnst = false
arrks_h.γ_cnst = false
errs_IDT,rate_IDT,h = AdjRRK.derv_test!(IDT_solver!,arrks_IDT,arrks_h,ts_IDT,rk4,h0,Nref)


#plots
labels = ["γ const" "proper lin"]
markers = [:circle :square]

errors = [errs_IDTγ0,errs_IDT]
plot(h,errors,
    title="linearization error (IDT)",
    xlabel="h",
    ylabel="error",
    label=labels,
    marker=markers,
    xaxis=:log,
    yaxis=:log)
display(plot!())

rates = [rate_IDTγ0,rate_IDT]
plot(rates,
    title="linearization conv rate (IDT)",
    xlabel="refinement index",
    ylabel="rate",
    label=labels,
    marker=markers,
    xaxis=:flip)
display(plot!())


##RRK
arrks_RRK.u0_lin = randn(2)

# assuming γ is cnst in linearization
arrks_RRK.γ_cnst = true
arrks_h.γ_cnst = true
errs_RRKγ0,rate_RRKγ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# with proper linearization, but const dt
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = true
arrks_h.γ_cnst = false
arrks_h.dt_cnst = true
errs_RRKdt0,rate_RRKdt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# with proper linearization
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = false
arrks_h.γ_cnst = false
arrks_h.dt_cnst = false
errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)


#plots
labels = ["γ const" "dt const" "proper lin"]
markers = [:circle :square :star5]

errors = [errs_RRKγ0,errs_RRKdt0,errs_RRK]
plot(h,errors,
    title="linearization error (RRK)",
    xlabel="h",
    ylabel="error",
    label=labels,
    marker=markers,
    xaxis=:log,
    yaxis=:log)
display(plot!())

rates = [rate_RRKγ0,rate_RRKdt0,rate_RRK]
plot(rates,
    title="linearization conv rate (RRK)",
    xlabel="refinement index",
    ylabel="rate",
    label=labels,
    marker=markers,
    xaxis=:flip)
display(plot!())
