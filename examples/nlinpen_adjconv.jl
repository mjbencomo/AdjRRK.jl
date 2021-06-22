# Testing out RRK method (linearized and adjoint code) on a small nonlinear ODE
# system (nlinear pendulum) from Ranocha 2019.

# Comparing the linearized and adjoint solutions from RRK algs when relaxation
# parameter γ is viewed as a constant or function of state fields during
# linearization. Also, looking at the effects of considering the last time step
# correction in the linearization.

using AdjRRK
using LinearAlgebra
using UnPack
using Plots
using LaTeXStrings
using Random
using MAT

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

#specifying 3rd order RK
s_rk = 3
A_rk = [0 0 0; 1 0 0; 1/4 1/4 0]
b_rk = [1/6; 1/6; 2/3]
rk3 = RK_struct(s_rk,b_rk,A_rk)

t0 = 0
T  = 2 #500
u0 = [1.5,1]
dt_small = 0.00001

Nref = 5
dt0 = 1
dt = zeros(Nref+1)
dt[1] = dt0
for n=1:Nref
    dt[n+1] = dt[n]/2
end

arrks_RK = AdjRRK_struct()
@pack! arrks_RK = f,df,η,∇η,Hη,u0
ts_RK = Time_struct()
@pack! ts_RK = t0,T

arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη,u0
ts_RRK = Time_struct()
@pack! ts_RRK = t0,T


file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK paper/figs/RRK_adj_conv.mat","w")
write(file,"dt",dt)
write(file,"T",T)
write(file,"dt_small",dt_small)


## computing "true" solution by RK with small dt
ts_RK.dt = dt_small

RK_solver!(arrks_RK,ts_RK,rk4)
arrks_RK.uT_adj = arrks_RK.u[:,end]
RK_solver!(arrks_RK,ts_RK,rk4;adj=true)

u_true = arrks_RK.u
z_true = arrks_RK.u_adj

norm_u = norm(u_true[:,end])
norm_z = norm(z_true[:,1])


## Convergence for RK2
rk = rk2

errs_fwd_RK2  = zeros(Nref+1)
errs_fwd_RRK2 = zeros(Nref+1)
errs_adj_RK2  = zeros(Nref+1)
errs_adj_RRK2 = zeros(Nref+1)
errs_adj_RRK2_g0  = zeros(Nref+1)
errs_adj_RRK2_dt0 = zeros(Nref+1)

for n=1:Nref+1
    ts_RK.dt  = dt[n]
    ts_RRK.dt = dt[n]

    #fwd RK
    RK_solver!(arrks_RK,ts_RK,rk)
    errs_fwd_RK2[n] = norm(arrks_RK.u[:,end]-u_true[:,end])/norm_u

    #adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    errs_adj_RK2[n] = norm(arrks_RK.u_adj[:,1]-z_true[:,1])/norm_z

    #fwd RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    errs_fwd_RRK2[n] = norm(arrks_RRK.u[:,end]-u_true[:,end])/norm_u

    #adj RRK
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK2[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK dt*-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK2_dt0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK γ-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK2_g0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z
end


rate_fwd_RK2  = zeros(Nref)
rate_fwd_RRK2 = zeros(Nref)
rate_adj_RK2  = zeros(Nref)
rate_adj_RRK2 = zeros(Nref)
rate_adj_RRK2_g0  = zeros(Nref)
rate_adj_RRK2_dt0 = zeros(Nref)

for n=1:Nref
    rate_fwd_RK2[n]  = log2(errs_fwd_RK2[n]) - log2(errs_fwd_RK2[n+1])
    rate_fwd_RRK2[n] = log2(errs_fwd_RRK2[n]) - log2(errs_fwd_RRK2[n+1])
    rate_adj_RK2[n]      = log2(errs_adj_RK2[n]) - log2(errs_adj_RK2[n+1])
    rate_adj_RRK2[n]     = log2(errs_adj_RRK2[n]) - log2(errs_adj_RRK2[n+1])
    rate_adj_RRK2_g0[n]  = log2(errs_adj_RRK2_g0[n]) - log2(errs_adj_RRK2_g0[n+1])
    rate_adj_RRK2_dt0[n] = log2(errs_adj_RRK2_dt0[n]) - log2(errs_adj_RRK2_dt0[n+1])
end

labels = ["RK" "RRK"]
markers = [:circle :square]
plot(dt,[errs_fwd_RK2,errs_fwd_RRK2],
    title="RK2 fwd error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK"]
markers = [:circle :square]
plot([rate_fwd_RK2,rate_fwd_RRK2],
    title="RK2 fwd rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot(dt,[errs_adj_RK2,errs_adj_RRK2,errs_adj_RRK2_g0,errs_adj_RRK2_dt0],
    title="RK2 adj error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot([rate_adj_RK2,rate_adj_RRK2,rate_adj_RRK2_g0,rate_adj_RRK2_dt0],
    title="RK2 adj rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


write(file,"errs_fwd_RK2",errs_fwd_RK2)
write(file,"errs_fwd_RRK2",errs_fwd_RRK2)
write(file,"errs_adj_RK2",errs_adj_RK2)
write(file,"errs_adj_RRK2",errs_adj_RRK2)
write(file,"errs_adj_RRK2_dt0",errs_adj_RRK2_dt0)
write(file,"errs_adj_RRK2_g0",errs_adj_RRK2_g0)


## Convergence for RK3
rk = rk3

errs_fwd_RK3  = zeros(Nref+1)
errs_fwd_RRK3 = zeros(Nref+1)
errs_adj_RK3  = zeros(Nref+1)
errs_adj_RRK3 = zeros(Nref+1)
errs_adj_RRK3_g0  = zeros(Nref+1)
errs_adj_RRK3_dt0 = zeros(Nref+1)

for n=1:Nref+1
    ts_RK.dt  = dt[n]
    ts_RRK.dt = dt[n]

    #fwd RK
    RK_solver!(arrks_RK,ts_RK,rk)
    errs_fwd_RK3[n] = norm(arrks_RK.u[:,end]-u_true[:,end])/norm_u

    #adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    errs_adj_RK3[n] = norm(arrks_RK.u_adj[:,1]-z_true[:,1])/norm_z

    #fwd RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    errs_fwd_RRK3[n] = norm(arrks_RRK.u[:,end]-u_true[:,end])/norm_u

    #adj RRK
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK3[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK dt*-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK3_dt0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK γ-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK3_g0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z
end


rate_fwd_RK3  = zeros(Nref)
rate_fwd_RRK3 = zeros(Nref)
rate_adj_RK3  = zeros(Nref)
rate_adj_RRK3 = zeros(Nref)
rate_adj_RRK3_g0  = zeros(Nref)
rate_adj_RRK3_dt0 = zeros(Nref)

for n=1:Nref
    rate_fwd_RK3[n]  = log2(errs_fwd_RK3[n]) - log2(errs_fwd_RK3[n+1])
    rate_fwd_RRK3[n] = log2(errs_fwd_RRK3[n]) - log2(errs_fwd_RRK3[n+1])
    rate_adj_RK3[n]      = log2(errs_adj_RK3[n]) - log2(errs_adj_RK3[n+1])
    rate_adj_RRK3[n]     = log2(errs_adj_RRK3[n]) - log2(errs_adj_RRK3[n+1])
    rate_adj_RRK3_g0[n]  = log2(errs_adj_RRK3_g0[n]) - log2(errs_adj_RRK3_g0[n+1])
    rate_adj_RRK3_dt0[n] = log2(errs_adj_RRK3_dt0[n]) - log2(errs_adj_RRK3_dt0[n+1])
end

labels = ["RK" "RRK"]
markers = [:circle :square]
plot(dt,[errs_fwd_RK3,errs_fwd_RRK3],
    title="RK3 fwd error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK"]
markers = [:circle :square]
plot([rate_fwd_RK3,rate_fwd_RRK3],
    title="RK3 fwd rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot(dt,[errs_adj_RK3,errs_adj_RRK3,errs_adj_RRK3_g0,errs_adj_RRK3_dt0],
    title="RK3 adj error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot([rate_adj_RK3,rate_adj_RRK3,rate_adj_RRK3_g0,rate_adj_RRK3_dt0],
    title="RK3 adj rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


write(file,"errs_fwd_RK3",errs_fwd_RK3)
write(file,"errs_fwd_RRK3",errs_fwd_RRK3)
write(file,"errs_adj_RK3",errs_adj_RK3)
write(file,"errs_adj_RRK3",errs_adj_RRK3)
write(file,"errs_adj_RRK3_dt0",errs_adj_RRK3_dt0)
write(file,"errs_adj_RRK3_g0",errs_adj_RRK3_g0)


## Convergence for RK4
rk = rk4

errs_fwd_RK4  = zeros(Nref+1)
errs_fwd_RRK4 = zeros(Nref+1)
errs_adj_RK4  = zeros(Nref+1)
errs_adj_RRK4 = zeros(Nref+1)
errs_adj_RRK4_g0  = zeros(Nref+1)
errs_adj_RRK4_dt0 = zeros(Nref+1)

for n=1:Nref+1
    ts_RK.dt  = dt[n]
    ts_RRK.dt = dt[n]

    #fwd RK
    RK_solver!(arrks_RK,ts_RK,rk)
    errs_fwd_RK4[n] = norm(arrks_RK.u[:,end]-u_true[:,end])/norm_u

    #adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    errs_adj_RK4[n] = norm(arrks_RK.u_adj[:,1]-z_true[:,1])/norm_z

    #fwd RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    errs_fwd_RRK4[n] = norm(arrks_RRK.u[:,end]-u_true[:,end])/norm_u

    #adj RRK
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK4[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK dt*-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK4_dt0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z

    #adj RRK γ-const
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    errs_adj_RRK4_g0[n] = norm(arrks_RRK.u_adj[:,1]-z_true[:,1])/norm_z
end


rate_fwd_RK4  = zeros(Nref)
rate_fwd_RRK4 = zeros(Nref)
rate_adj_RK4  = zeros(Nref)
rate_adj_RRK4 = zeros(Nref)
rate_adj_RRK4_g0  = zeros(Nref)
rate_adj_RRK4_dt0 = zeros(Nref)

for n=1:Nref
    rate_fwd_RK4[n]  = log2(errs_fwd_RK4[n]) - log2(errs_fwd_RK4[n+1])
    rate_fwd_RRK4[n] = log2(errs_fwd_RRK4[n]) - log2(errs_fwd_RRK4[n+1])
    rate_adj_RK4[n]      = log2(errs_adj_RK4[n]) - log2(errs_adj_RK4[n+1])
    rate_adj_RRK4[n]     = log2(errs_adj_RRK4[n]) - log2(errs_adj_RRK4[n+1])
    rate_adj_RRK4_g0[n]  = log2(errs_adj_RRK4_g0[n]) - log2(errs_adj_RRK4_g0[n+1])
    rate_adj_RRK4_dt0[n] = log2(errs_adj_RRK4_dt0[n]) - log2(errs_adj_RRK4_dt0[n+1])
end

labels = ["RK" "RRK"]
markers = [:circle :square]
plot(dt,[errs_fwd_RK4,errs_fwd_RRK4],
    title="RK4 fwd error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK"]
markers = [:circle :square]
plot([rate_fwd_RK4,rate_fwd_RRK4],
    title="RK4 fwd rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot(dt,[errs_adj_RK4,errs_adj_RRK4,errs_adj_RRK4_g0,errs_adj_RRK4_dt0],
    title="RK4 adj error",
    xaxis=:log,
    yaxis=:log,
    label=labels,
    marker=markers,
    legend=:topleft,
    xlabel="dt")
display(plot!())

labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
markers = [:circle :square :cross :star]
plot([rate_adj_RK4,rate_adj_RRK4,rate_adj_RRK4_g0,rate_adj_RRK4_dt0],
    title="RK4 adj rates",
    label=labels,
    marker=markers,
    legend=:bottomright)
display(plot!())


write(file,"errs_fwd_RK4",errs_fwd_RK4)
write(file,"errs_fwd_RRK4",errs_fwd_RRK4)
write(file,"errs_adj_RK4",errs_adj_RK4)
write(file,"errs_adj_RRK4",errs_adj_RRK4)
write(file,"errs_adj_RRK4_dt0",errs_adj_RRK4_dt0)
write(file,"errs_adj_RRK4_g0",errs_adj_RRK4_g0)


close(file)
