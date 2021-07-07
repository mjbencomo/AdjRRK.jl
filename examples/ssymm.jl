# This code runs some numerical experiments and generates .mat files for adjoint
# RRK paper. In particular, we run several experiments for a random skew symm
# initial value problem and study the convergence of the time-reverse errors in
# which we solve the forward problem and then the adjoint problem and see if we
# get back the initial condition.

using AdjRRK
using LinearAlgebra
using Plots
using UnPack
using LaTeXStrings
using Random
using MAT

N = 10
Random.seed!(1234)
A = randn(N,N)
A = A - transpose(A)

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


function time_symm_test!(solver!,arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct,dt0,Nref,true_sol)
    dt = zeros(Nref+1)
    dt[1] = dt0
    for n=1:Nref
        dt[n+1] = dt[n]/2
    end

    norm_true_sol = norm(true_sol)

    errs = zeros(Nref+1)
    for n=1:Nref+1
        ts.dt  = dt[n]
        solver!(arrks,ts,rk)
        arrks.uT_adj = arrks.u[:,end]
        solver!(arrks,ts,rk;adj=true)

        errs[n] = norm(arrks.u_adj[:,1]-true_sol)/norm_true_sol
    end

    rate = zeros(Nref)
    for n=1:Nref
        rate[n]  = log2(errs[n]) - log2(errs[n+1])
    end

    return errs, rate, dt
end


# Flags
write_mat = true  # Want to output .mat file?
make_plot = true # Want to output plots?

if write_mat
    file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK paper/figs/RRK_ssymm.mat","w")
end

# Initial time and initial condition
t0  = 0
T   = norm(A)*10
Random.seed!(12345)
u0 = randn(N)

# Initializing RRK and time axis structs
arrks = AdjRRK_struct()
@pack! arrks = η,∇η,Hη,df
arrks.f = A
@pack! arrks = u0


ts = Time_struct()
@pack! ts = t0,T

# Refinement info
dt0 = 1/norm(A)
dt0 = dt0/2
Nref = 6


# DIRK3 case ------------------------------------------------------------------#
rk = dirk3

# RK
errs_RK,rate_RK,dt = time_symm_test!(DIRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK
arrks.dt_cnst = false
arrks.γ_cnst  = false
errs_RRK,rate_RRK,dt = time_symm_test!(DIRRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (dt*-const)
arrks.dt_cnst = true
arrks.γ_cnst  = false
errs_RRK_dt0,rate_RRK_dt0,dt = time_symm_test!(DIRRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (γ-const)
arrks.dt_cnst = true
arrks.γ_cnst  = true
errs_RRK_g0,rate_RRK_g0,dt = time_symm_test!(DIRRK_solver!,arrks,ts,rk,dt0,Nref,u0)

if make_plot
    labels = ["DIRK3" "DIRRK3" "DIRRK3 (Δt*-const)" "DIRRK3 (γ-const)"]
    markers = [:circle :square :square :star]
    linestyles = [:dash :dash :dash :dash]

    plot(dt,[errs_RK,errs_RRK,errs_RRK_dt0,errs_RRK_g0],
        xaxis=:log,
        yaxis=:log,
        label=labels,
        legend=:topleft,
        marker=markers,
        linestyle=linestyles,
        xlabel=L"\Delta t",
        ylabel="Error at t=$T")
    display(plot!())

    plot([rate_RK,rate_RRK,rate_RRK_dt0,rate_RRK_g0],
        label=labels,
        marker=markers,
        linestyle=linestyles,
        xaxis=:flip,
        title="Adjoint convergence slopes/rates",
        ylabel="Rates/slopes",
        xlabel="refinement index")
    display(plot!())
end

if write_mat
    write(file,"errs_DIRK3",errs_RK)
    write(file,"errs_DIRRK3",errs_RRK)
    write(file,"errs_DIRRK3_dt0",errs_RRK_dt0)
    write(file,"errs_DIRRK3_g0",errs_RRK_g0)

    write(file,"rate_DIRK3",rate_RK)
    write(file,"rate_DIRRK3",rate_RRK)
    write(file,"rate_DIRRK3_dt0",rate_RRK_dt0)
    write(file,"rate_DIRRK3_g0",rate_RRK_g0)
end


# Need to add this for explicit RK
@pack! arrks = f

# RK4 case --------------------------------------------------------------------#
rk = rk4

# RK
errs_RK,rate_RK,dt = time_symm_test!(RK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK
arrks.dt_cnst = false
arrks.γ_cnst  = false
errs_RRK,rate_RRK,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (dt*-const)
arrks.dt_cnst = true
arrks.γ_cnst  = false
errs_RRK_dt0,rate_RRK_dt0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (γ-const)
arrks.dt_cnst = true
arrks.γ_cnst  = true
errs_RRK_g0,rate_RRK_g0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

if make_plot
    labels = ["RK4" "RRK4" "RRK4 (Δt*-const)" "RRK4 (γ-const)"]
    markers = [:circle :square :square :star]
    linestyles = [:dash :dash :dash :dash]

    plot(dt,[errs_RK,errs_RRK,errs_RRK_dt0,errs_RRK_g0],
        xaxis=:log,
        yaxis=:log,
        label=labels,
        legend=:topleft,
        marker=markers,
        linestyle=linestyles,
        xlabel=L"\Delta t",
        ylabel="Error at t=$T")
    display(plot!())

    plot([rate_RK,rate_RRK,rate_RRK_dt0,rate_RRK_g0],
        label=labels,
        marker=markers,
        linestyle=linestyles,
        xaxis=:flip,
        title="Adjoint convergence slopes/rates",
        ylabel="Rates/slopes",
        xlabel="refinement index")
    display(plot!())
end

if write_mat
    write(file,"errs_RK4",errs_RK)
    write(file,"errs_RRK4",errs_RRK)
    write(file,"errs_RRK4_dt0",errs_RRK_dt0)
    write(file,"errs_RRK4_g0",errs_RRK_g0)

    write(file,"rate_RK4",rate_RK)
    write(file,"rate_RRK4",rate_RRK)
    write(file,"rate_RRK4_dt0",rate_RRK_dt0)
    write(file,"rate_RRK4_g0",rate_RRK_g0)
end


# RK3 case --------------------------------------------------------------------#
rk = rk3

# RK
errs_RK,rate_RK,dt = time_symm_test!(RK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK
arrks.dt_cnst = false
arrks.γ_cnst  = false
errs_RRK,rate_RRK,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (dt*-const)
arrks.dt_cnst = true
arrks.γ_cnst  = false
errs_RRK_dt0,rate_RRK_dt0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (γ-const)
arrks.dt_cnst = true
arrks.γ_cnst  = true
errs_RRK_g0,rate_RRK_g0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

if make_plot
    labels = ["RK3" "RRK3" "RRK3 (Δt*-const)" "RRK3 (γ-const)"]
    markers = [:circle :square :square :star]
    linestyles = [:dash :dash :dash :dash]

    plot(dt,[errs_RK,errs_RRK,errs_RRK_dt0,errs_RRK_g0],
        xaxis=:log,
        yaxis=:log,
        label=labels,
        legend=:topleft,
        marker=markers,
        linestyle=linestyles,
        xlabel=L"\Delta t",
        ylabel="Error at t=$T")
    display(plot!())

    plot([rate_RK,rate_RRK,rate_RRK_dt0,rate_RRK_g0],
        label=labels,
        marker=markers,
        linestyle=linestyles,
        xaxis=:flip,
        title="Adjoint convergence slopes/rates",
        ylabel="Rates/slopes",
        xlabel="refinement index")
    display(plot!())
end

if write_mat
    write(file,"errs_RK3",errs_RK)
    write(file,"errs_RRK3",errs_RRK)
    write(file,"errs_RRK3_dt0",errs_RRK_dt0)
    write(file,"errs_RRK3_g0",errs_RRK_g0)

    write(file,"rate_RK3",rate_RK)
    write(file,"rate_RRK3",rate_RRK)
    write(file,"rate_RRK3_dt0",rate_RRK_dt0)
    write(file,"rate_RRK3_g0",rate_RRK_g0)
end


# RK2 case --------------------------------------------------------------------#
rk = rk2

# RK
errs_RK,rate_RK,dt = time_symm_test!(RK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK
arrks.dt_cnst = false
arrks.γ_cnst  = false
errs_RRK,rate_RRK,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (dt*-const)
arrks.dt_cnst = true
arrks.γ_cnst  = false
errs_RRK_dt0,rate_RRK_dt0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

# RRK (γ-const)
arrks.dt_cnst = true
arrks.γ_cnst  = true
errs_RRK_g0,rate_RRK_g0,dt = time_symm_test!(RRK_solver!,arrks,ts,rk,dt0,Nref,u0)

if make_plot
    labels = ["RK2" "RRK2" "RRK2 (Δt*-const)" "RRK2 (γ-const)"]
    markers = [:circle :square :square :star]
    linestyles = [:dash :dash :dash :dash]

    plot(dt,[errs_RK,errs_RRK,errs_RRK_dt0,errs_RRK_g0],
        xaxis=:log,
        yaxis=:log,
        label=labels,
        legend=:topleft,
        marker=markers,
        linestyle=linestyles,
        xlabel=L"\Delta t",
        ylabel="Error at t=$T")
    display(plot!())

    plot([rate_RK,rate_RRK,rate_RRK_dt0,rate_RRK_g0],
        label=labels,
        marker=markers,
        linestyle=linestyles,
        xaxis=:flip,
        title="Adjoint convergence slopes/rates",
        ylabel="Rates/slopes",
        xlabel="refinement index")
    display(plot!())
end

if write_mat
    write(file,"errs_RK2",errs_RK)
    write(file,"errs_RRK2",errs_RRK)
    write(file,"errs_RRK2_dt0",errs_RRK_dt0)
    write(file,"errs_RRK2_g0",errs_RRK_g0)

    write(file,"rate_RK2",rate_RK)
    write(file,"rate_RRK2",rate_RRK)
    write(file,"rate_RRK2_dt0",rate_RRK_dt0)
    write(file,"rate_RRK2_g0",rate_RRK_g0)
end

if write_mat
    write(file,"dt",dt)
    close(file)
end
