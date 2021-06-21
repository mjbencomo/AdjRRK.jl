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
# using Plots.Measures

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
dt = 0.9
u0 = [1.5,1]

#Plotting params
# fntsm = Plots.font("sans-serif", pointsize=round(10))
# fntlg = Plots.font("sans-serif", pointsize=round(15))
# default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)


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

plot(title="RK solution",
    u1_RK,u2_RK,
    xlabel=L"y_1",
    ylabel=L"y_2",
    legend=false)
    # xlims=(-2,2),
    # xticks=[-2,-1,0,1,2],
    # yticks=[-2,-1,0,1,2],
    # linewidth=0.1)
# plot!(size(1000,1000))
#plot!(left_margin=5mm, bottom_margin=5mm, right_margin=5mm)
display(plot!())
# savefig("RK_soln.pdf")

plot(title="RK entropy production",
    ts_RK.t,arrks_RK.Δη,
    xlabel=L"t",
    ylabel=L"\eta(t)-\eta(0)",
    legend=false)
display(plot!())
# savefig("RK_entropy.pdf")


#running lin RK
arrks_RK.u0_lin = u0
RK_solver!(arrks_RK,ts_RK,rk4;lin=true)
w1_RK = arrks_RK.u_lin[1,:]
w2_RK = arrks_RK.u_lin[2,:]

plot(title="linearized RK solution",
    w1_RK,w2_RK,
    xlabel=L"\delta_1",
    ylabel=L"\delta_2",
    legend=false)
display(plot!())
# savefig("RK_soln_lin.pdf")


#running adj RK
arrks_RK.uT_adj = arrks_RK.u[:,end]
RK_solver!(arrks_RK,ts_RK,rk4;adj=true)
z1_RK = arrks_RK.u_adj[1,:]
z2_RK = arrks_RK.u_adj[2,:]

plot(title="adjoint RK solution",
    z1_RK,z2_RK,
    xlabel=L"\lambda_1",
    ylabel=L"\lambda_2",
    legend=false)
display(plot!())
# savefig("RK_soln_adj.pdf")


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

plot(title="RRK solution",
    u1_RRK,u2_RRK,
    xlabel=L"y_1",
    ylabel=L"y_2",
    legend=false)
display(plot!())
# savefig("RRK_soln.pdf")

plot(title="RRK entropy production",
    ts_RRK.t,arrks_RRK.Δη,
    xlabel=L"t",
    ylabel=L"\eta(t)-\eta(0)",
    legend=false)
    # yticks=[0,5e-15,1e-14,1.5e-14]
    #)
display(plot!())
# savefig("RRK_entropy.pdf")

plot(title="Relaxation parameter",
    ts_RRK.t,arrks_RRK.γ,
    xlabel=L"t",
    ylabel=L"\gamma",
    legend=false)
    # yticks=[0,5e-15,1e-14,1.5e-14]
    #)
display(plot!())

#running lin RRK (γ0)
arrks_RRK.γ_cnst = true
arrks_RRK.u0_lin = u0
RRK_solver!(arrks_RRK,ts_RRK,rk4;lin=true)
w1_RRKγ0 = arrks_RRK.u_lin[1,:]
w2_RRKγ0 = arrks_RRK.u_lin[2,:]

plot(title="linearized RRK (γ-const) solution",
    w1_RRKγ0,w2_RRKγ0,
    xlabel=L"\delta_1",
    ylabel=L"\delta_2",
    legend=false)
display(plot!())
# savefig("RRK0_soln_lin.pdf")

#running lin RRK
arrks_RRK.γ_cnst = false
arrks_RRK.u0_lin = u0
RRK_solver!(arrks_RRK,ts_RRK,rk4;lin=true)
w1_RRK = arrks_RRK.u_lin[1,:]
w2_RRK = arrks_RRK.u_lin[2,:]

plot(title="linearized RRK solution",
    w1_RRK,w2_RRK,
    xlabel=L"\delta_1",
    ylabel=L"\delta_2",
    legend=false)
display(plot!())
# savefig("RRK_soln_lin.pdf")


#running adj RRK (γ-constant)
arrks_RRK.γ_cnst = true
arrks_RRK.dt_cnst = true
arrks_RRK.uT_adj = arrks_RRK.u[:,end]
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
z1_RRKγ0 = arrks_RRK.u_adj[1,:]
z2_RRKγ0 = arrks_RRK.u_adj[2,:]

plot(title="adjoint RRK (γ-const) solution",
    z1_RRKγ0,z2_RRKγ0,
    xlabel=L"\lambda_1",
    ylabel=L"\lambda_2",
    legend=false)
display(plot!())
# savefig("RRK0_soln_adj.pdf")


#running adj RRK (Δt*-constant)
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = true
arrks_RRK.uT_adj = arrks_RRK.u[:,end]
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
z1_RRKdt0 = arrks_RRK.u_adj[1,:]
z2_RRKdt0 = arrks_RRK.u_adj[2,:]

plot(title="adjoint RRK (Δt*-const) solution",
    z1_RRKdt0,z2_RRKdt0,
    xlabel=L"\lambda_1",
    ylabel=L"\lambda_2",
    legend=false)
display(plot!())


#running adj RRK
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = false
arrks_RRK.uT_adj = arrks_RRK.u[:,end]
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
z1_RRK = arrks_RRK.u_adj[1,:]
z2_RRK = arrks_RRK.u_adj[2,:]

plot(title="adjoint RRK solution",
    z1_RRK,z2_RRK,
    xlabel=L"\lambda_1",
    ylabel=L"\lambda_2",
    legend=false)
display(plot!())
# savefig("RRK_soln_adj.pdf")


##derivative tests

Nref = 5
h0 = 2^(-12)

arrks_h = AdjRRK_struct()
@pack! arrks_h = f,df,η,∇η,Hη


# ##IDT
# arrks_IDT.u0_lin = randn(2)
#
# # assuming γ is cnst in linearization
# arrks_IDT.γ_cnst = true
# arrks_h.γ_cnst = true
# errs_IDTγ0,rate_IDTγ0,h = AdjRRK.derv_test!(IDT_solver!,arrks_IDT,arrks_h,ts_IDT,rk4,h0,Nref)
#
# # with proper linearization
# arrks_IDT.γ_cnst = false
# arrks_h.γ_cnst = false
# errs_IDT,rate_IDT,h = AdjRRK.derv_test!(IDT_solver!,arrks_IDT,arrks_h,ts_IDT,rk4,h0,Nref)
#
#
# #plots
# labels = ["γ const" "proper lin"]
# markers = [:circle :square]
#
# errors = [errs_IDTγ0,errs_IDT]
# plot(h,errors,
#     title="linearization error (IDT)",
#     xlabel="h",
#     ylabel="error",
#     label=labels,
#     marker=markers,
#     xaxis=:log,
#     yaxis=:log)
# display(plot!())
# savefig("linerr_IDT.png")
#
# rates = [rate_IDTγ0,rate_IDT]
# plot(rates,
#     title="linearization conv rate (IDT)",
#     xlabel="refinement index",
#     ylabel="rate",
#     label=labels,
#     marker=markers,
#     xaxis=:flip)
# display(plot!())


##RRK
Random.seed!(123)
arrks_RRK.u0_lin = randn(2)

# assuming γ is cnst in linearization
arrks_RRK.γ_cnst = true
arrks_RRK.dt_cnst = true
# arrks_h.γ_cnst = true
# arrks_h.dt_cnst = true
errs_RRKγ0,rate_RRKγ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# with proper linearization, but const dt
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = true
# arrks_h.γ_cnst = false
# arrks_h.dt_cnst = true
errs_RRKdt0,rate_RRKdt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# with proper linearization
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = false
# arrks_h.γ_cnst = false
# arrks_h.dt_cnst = false
errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)


#plots
labels = ["γ constant" "Δt* constant" "proper lin."]
markers = [:circle :square :star5]

errors = [errs_RRKγ0,errs_RRKdt0,errs_RRK]
plot(title="linearization error",
    h,errors,
    xlabel=L"h",
    ylabel="error",
    label=labels,
    marker=markers,
    xaxis=:log,
    yaxis=:log)
display(plot!())
# savefig("linerr_RRK.pdf")

# rates = [rate_RRKγ0,rate_RRKdt0,rate_RRK]
# plot(rates,
#     title="linearization conv rate (RRK)",
#     xlabel="refinement index",
#     ylabel="rate",
#     label=labels,
#     marker=markers,
#     xaxis=:flip)
# display(plot!())
