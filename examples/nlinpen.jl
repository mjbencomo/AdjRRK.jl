# This code runs some numerical experiments and generates .mat files for adjoint
# RRK paper. In particular, we run several experiments for a simple nonlinear
# pendulum problem (Ranocha 2019):
#   1. Derivative tests for verifying linearization of RRK and quantifying
#      errors from improper linearization.
#   2. Convergence experiment of discrete adjoint.
#   3. Accuracy experiment of discrete adjoint.

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


function adj_conv_test!(solver!,arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct,dt0,Nref,true_sol,fwd_flag)
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

        if !fwd_flag
            arrks.uT_adj = arrks.u[:,end]
            solver!(arrks,ts,rk;adj=true)
        end

        if fwd_flag
            errs[n] = norm(arrks.u[:,end]-true_sol)/norm_true_sol
        else
            errs[n] = norm(arrks.u_adj[:,1]-true_sol)/norm_true_sol
        end
    end

    rate = zeros(Nref)
    for n=1:Nref
        rate[n]  = log2(errs[n]) - log2(errs[n+1])
    end

    return errs, rate, dt
end


# Flags
write_mat = true  # Want to output .mat file?
make_plot = false # Want to output plots?

run_derv_test  = false # Want to run derivative tests?
run_adj_conv   = false # Want to run adjoint convergence tests?
run_adj_growth = true # Want to run adjoint growth tests?


# Initial time and initial condition
t0 = 0
u0 = [1.5,1]

# Initializing RRK and time axis structs
arrks_RK = AdjRRK_struct()
@pack! arrks_RK = f,df,η,∇η,Hη
@pack! arrks_RK = u0
arrks_RK.return_time = true
arrks_RK.return_Δη = true

arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη
@pack! arrks_RRK = u0
arrks_RRK.return_time = true
arrks_RRK.return_Δη = true

ts_RK = Time_struct()
@pack! ts_RK = t0
ts_RRK = Time_struct()
@pack! ts_RRK = t0


#==============================================================================#
# 1. Derivative test: checking that a first order FD approximation to
#    directional derivative converges appropriately.

if run_derv_test
    print("\nRunning derivative tests .... ")

    if write_mat
        file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK_paper/figs/RRK_nlpen_derv.mat","w")
    end

    ts_RK.T  = 200
    ts_RRK.T = 200

    dt_RK4 = 0.1#0.9
    dt_RK3 = 0.1#0.65
    dt_RK2 = 0.1#0.4

    Nref = 12
    h0 = 2^(-12)

    arrks_h = AdjRRK_struct()
    @pack! arrks_h = f,df,η,∇η,Hη

    Random.seed!(123)
    arrks_RRK.u0_lin = randn(2)


    # RK4 case ----------------------------------------------------------------#
    rk = rk4
    ts_RRK.dt = dt_RK4

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["proper lin." "Δt*-const" "γ-const"]
        markers = [:circle :square :star5]

        errors = [errs_RRK,errs_RRK_dt0,errs_RRK_γ0]
        plot(title="FD error (RRK4)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK,rate_RRK_dt0,rate_RRK_γ0]
        # plot(rates,
        #     title="FD convergence rate (RRK4)",
        #     xlabel="refinement index",
        #     ylabel="rate",
        #     label=labels,
        #     marker=markers,
        #     xaxis=:flip)
        # display(plot!())
    end

    if write_mat
        write(file,"errs_RRK4",errs_RRK)
        write(file,"errs_RRK4_dt0",errs_RRK_dt0)
        write(file,"errs_RRK4_g0",errs_RRK_γ0)
    end


    # RK3 case ----------------------------------------------------------------#
    rk = rk3
    ts_RRK.dt = dt_RK3

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["proper lin." "Δt* constant" "γ constant"]
        markers = [:circle :square :star5]

        errors = [errs_RRK,errs_RRK_dt0,errs_RRK_γ0]
        plot(title="FD error (RK3)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK,rate_RRK_dt0,rate_RRK_γ0]
        # plot(rates,
        #     title="FD convergence rate (RRK3)",
        #     xlabel="refinement index",
        #     ylabel="rate",
        #     label=labels,
        #     marker=markers,
        #     xaxis=:flip)
        # display(plot!())
    end

    if write_mat
        write(file,"errs_RRK3",errs_RRK)
        write(file,"errs_RRK3_dt0",errs_RRK_dt0)
        write(file,"errs_RRK3_g0",errs_RRK_γ0)
    end


    # RK2 case ----------------------------------------------------------------#
    rk = rk2
    ts_RRK.dt = dt_RK2

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["proper lin." "Δt* constant" "γ constant"]
        markers = [:circle :square :star5]

        errors = [errs_RRK,errs_RRK_dt0,errs_RRK_γ0]
        plot(title="FD error (RK2)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK,rate_RRK_dt0,rate_RRK_γ0]
        # plot(rates,
        #     title="FD convergence rate (RRK)",
        #     xlabel="refinement index",
        #     ylabel="rate",
        #     label=labels,
        #     marker=markers,
        #     xaxis=:flip)
        # display(plot!())
    end

    if write_mat
        write(file,"errs_RRK2",errs_RRK)
        write(file,"errs_RRK2_dt0",errs_RRK_dt0)
        write(file,"errs_RRK2_g0",errs_RRK_γ0)
    end


    if write_mat
        write(file,"h_derv_test",h)
        close(file)
    end
    println("End of derivative tests.")
end


#==============================================================================#
# 2. Convergence of discrete adjoint.

if run_adj_conv
    print("\nRunning adjoint convergence tests .... ")
    if write_mat
        file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK_paper/figs/RRK_nlpen_conv.mat","w")
    end

    ts_RK.T = 2
    ts_RRK.T = 2
    dt_small = 0.00001

    dt0 = 1
    Nref = 5

    if write_mat
        write(file,"T_adj_conv",ts_RK.T)
        write(file,"dt_small",dt_small)
    end

    # computing "true" solution by RK with small dt
    ts_RK.dt = dt_small
    RK_solver!(arrks_RK,ts_RK,rk4)
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk4;adj=true)
    uT_true = arrks_RK.u[:,end]
    z0_true = arrks_RK.u_adj[:,1]


    # RK4 case ----------------------------------------------------------------#
    rk = rk4

    # fwd RK case
    fwd_flag = true
    cerrs_fwd_RK,crate_fwd_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,uT_true,fwd_flag)

    # fwd RRK case
    fwd_flag = true
    cerrs_fwd_RRK,crate_fwd_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,uT_true,fwd_flag)

    # adj RK case
    fwd_flag = false
    cerrs_adj_RK,crate_adj_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case
    fwd_flag = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK,crate_adj_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (dt*-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK_dt0,crate_adj_RRK_dt0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (γ-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    cerrs_adj_RRK_γ0,crate_adj_RRK_γ0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    if make_plot
        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot(dt,[cerrs_fwd_RK, cerrs_fwd_RRK],
            title="Forward error (RK4)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot([crate_fwd_RK, crate_fwd_RRK],
            title="Forward rates (RK4)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())

        labels = ["RK" "RRK" "RRK (Δt*-const)" "RRK (γ-const)"]
        markers = [:circle :square :cross :star]
        plot(dt,[cerrs_adj_RK, cerrs_adj_RRK, cerrs_adj_RRK_dt0, cerrs_adj_RRK_γ0],
            title="Adjoint error (RK4)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK" "RRK (Δt*-const)" "RRK (γ-const)"]
        markers = [:circle :square :cross :star]
        plot([crate_adj_RK, crate_adj_RRK, crate_adj_RRK_dt0, crate_adj_RRK_γ0],
            title="Adjoint rates (RK4)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())
    end

    if write_mat
        write(file,"cerrs_fwd_RK4",cerrs_fwd_RK)
        write(file,"cerrs_fwd_RRK4",cerrs_fwd_RRK)
        write(file,"cerrs_adj_RK4",cerrs_adj_RK)
        write(file,"cerrs_adj_RRK4",cerrs_adj_RRK)
        write(file,"cerrs_adj_RRK4_dt0",cerrs_adj_RRK_dt0)
        write(file,"cerrs_adj_RRK4_g0",cerrs_adj_RRK_γ0)
    end


    # RK3 ---------------------------------------------------------------------#
    rk = rk3

    # fwd RK case
    fwd_flag = true
    cerrs_fwd_RK,crate_fwd_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,uT_true,fwd_flag)

    # fwd RRK case
    fwd_flag = true
    cerrs_fwd_RRK,crate_fwd_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,uT_true,fwd_flag)

    # adj RK case
    fwd_flag = false
    cerrs_adj_RK,crate_adj_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case
    fwd_flag = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK,crate_adj_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (dt*-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK_dt0,crate_adj_RRK_dt0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (γ-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    cerrs_adj_RRK_γ0,crate_adj_RRK_γ0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    if make_plot
        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot(dt,[cerrs_fwd_RK, cerrs_fwd_RRK],
            title="Forward error (RK3)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot([crate_fwd_RK, crate_fwd_RRK],
            title="Forward rates (RK3)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())

        labels = ["RK" "RRK" "RRK (Δt*-const)" "RRK (γ-const)"]
        markers = [:circle :square :cross :star]
        plot(dt,[cerrs_adj_RK, cerrs_adj_RRK, cerrs_adj_RRK_dt0, cerrs_adj_RRK_γ0],
            title="Adjoint error (RK3)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK" "RRK (Δt*-const)" "RRK (γ-const)" ]
        markers = [:circle :square :cross :star]
        plot([crate_adj_RK, crate_adj_RRK, crate_adj_RRK_dt0, crate_adj_RRK_γ0],
            title="Adjoint rates (RK3)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())
    end

    if write_mat
        write(file,"cerrs_fwd_RK3",cerrs_fwd_RK)
        write(file,"cerrs_fwd_RRK3",cerrs_fwd_RRK)
        write(file,"cerrs_adj_RK3",cerrs_adj_RK)
        write(file,"cerrs_adj_RRK3",cerrs_adj_RRK)
        write(file,"cerrs_adj_RRK3_dt0",cerrs_adj_RRK_dt0)
        write(file,"cerrs_adj_RRK3_g0",cerrs_adj_RRK_γ0)
    end


    # RK2 ---------------------------------------------------------------------#
    rk = rk2

    # fwd RK case
    fwd_flag = true
    cerrs_fwd_RK,crate_fwd_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,uT_true,fwd_flag)

    # fwd RRK case
    fwd_flag = true
    cerrs_fwd_RRK,crate_fwd_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,uT_true,fwd_flag)

    # adj RK case
    fwd_flag = false
    cerrs_adj_RK,crate_adj_RK,dt = adj_conv_test!(RK_solver!,arrks_RK,ts_RK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case
    fwd_flag = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK,crate_adj_RRK,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (dt*-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = false
    cerrs_adj_RRK_dt0,crate_adj_RRK_dt0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    # adj RRK case (γ-const)
    fwd_flag = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.γ_cnst  = true
    cerrs_adj_RRK_γ0,crate_adj_RRK_γ0,dt = adj_conv_test!(RRK_solver!,arrks_RRK,ts_RRK,rk,dt0,Nref,z0_true,fwd_flag)

    if make_plot
        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot(dt,[cerrs_fwd_RK, cerrs_fwd_RRK],
            title="Forward error (RK2)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK"]
        markers = [:circle :square]
        plot([crate_fwd_RK, crate_fwd_RRK],
            title="Forward rates (RK2)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())

        labels = ["RK" "RRK" "RRK (γ-const)" "RRK (Δt*-const)"]
        markers = [:circle :square :cross :star]
        plot(dt,[cerrs_adj_RK, cerrs_adj_RRK, cerrs_adj_RRK_dt0, cerrs_adj_RRK_γ0],
            title="Adjoint error (RK2)",
            xaxis=:log,
            yaxis=:log,
            label=labels,
            marker=markers,
            legend=:topleft,
            xlabel="dt")
        display(plot!())

        labels = ["RK" "RRK" "RRK (Δt*-const)" "RRK (γ-const)"]
        markers = [:circle :square :cross :star]
        plot([crate_adj_RK, crate_adj_RRK, crate_adj_RRK_dt0, crate_adj_RRK_γ0],
            title="Adjoint rates (RK2)",
            label=labels,
            marker=markers,
            legend=:bottomright)
        display(plot!())
    end

    if write_mat
        write(file,"cerrs_fwd_RK2",cerrs_fwd_RK)
        write(file,"cerrs_fwd_RRK2",cerrs_fwd_RRK)
        write(file,"cerrs_adj_RK2",cerrs_adj_RK)
        write(file,"cerrs_adj_RRK2",cerrs_adj_RRK)
        write(file,"cerrs_adj_RRK2_g0",cerrs_adj_RRK_γ0)
        write(file,"cerrs_adj_RRK2_dt0",cerrs_adj_RRK_dt0)
    end


    if write_mat
        write(file,"dt_adj_conv",dt)
        close(file)
    end

    println("End of adjoint convergence tests.")
end


#==============================================================================#
# 3. Accuracy of discrete adjoint experiment.

if run_adj_growth
    print("\nRunning adjoint growth tests ... ")
    if write_mat
        file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK_paper/figs/RRK_nlpen_adj.mat","w")
    end

    ts_RK.T  = 200
    ts_RRK.T = 200

    dt_RK4 = 0.1
    dt_RK3 = 0.1
    dt_RK2 = 0.1


    # RK4 case ----------------------------------------------------------------#
    rk = rk4
    ts_RK.dt = dt_RK4
    ts_RRK.dt = dt_RK4

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)

    if make_plot
        plot(title="RK4 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK4 solution",
            arrks_RK.u[1,:],
            arrks_RK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK4 adjoint solution",
            arrks_RK.u_adj[1,:],
            arrks_RK.u_adj[2,:],
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK4",dt_RK4)
        write(file,"t_RK4",ts_RK.t)
        write(file,"y1_fwd_RK4",arrks_RK.u[1,:])
        write(file,"y2_fwd_RK4",arrks_RK.u[2,:])
        write(file,"y1_adj_RK4",arrks_RK.u_adj[1,:])
        write(file,"y2_adj_RK4",arrks_RK.u_adj[2,:])
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1 = arrks_RRK.u_adj[1,:]
    z2 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_dt0 = arrks_RRK.u_adj[1,:]
    z2_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_γ0 = arrks_RRK.u_adj[1,:]
    z2_γ0 = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK4 solution",
            arrks_RRK.u[1,:],
            arrks_RRK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 entropy production",
            ts_RRK.t,arrks_RRK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RRK4 relaxation parameter",
            ts_RRK.t,arrks_RRK.γ,
            xlabel=L"t",
            ylabel=L"\gamma",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution",
            z1,z2,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (Δt*-const)",
            z1_dt0,z2_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (γ-const)",
            z1_γ0,z2_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK4",ts_RRK.t)
        write(file,"y1_fwd_RRK4",arrks_RRK.u[1,:])
        write(file,"y2_fwd_RRK4",arrks_RRK.u[2,:])
        write(file,"y1_adj_RRK4",z1)
        write(file,"y2_adj_RRK4",z2)
        write(file,"y1_adj_RRK4_dt0",z1_dt0)
        write(file,"y2_adj_RRK4_dt0",z2_dt0)
        write(file,"y1_adj_RRK4_g0",z1_γ0)
        write(file,"y2_adj_RRK4_g0",z2_γ0)
    end


    # RK3 case ----------------------------------------------------------------#
    rk = rk3
    ts_RK.dt = dt_RK3
    ts_RRK.dt = dt_RK3

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)

    if make_plot
        plot(title="RK3 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK3 solution",
            arrks_RK.u[1,:],
            arrks_RK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK3 adjoint solution",
            arrks_RK.u_adj[1,:],
            arrks_RK.u_adj[2,:],
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK3",dt_RK3)
        write(file,"t_RK3",ts_RK.t)
        write(file,"y1_fwd_RK3",arrks_RK.u[1,:])
        write(file,"y2_fwd_RK3",arrks_RK.u[2,:])
        write(file,"y1_adj_RK3",arrks_RK.u_adj[1,:])
        write(file,"y2_adj_RK3",arrks_RK.u_adj[2,:])
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1 = arrks_RRK.u_adj[1,:]
    z2 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_dt0 = arrks_RRK.u_adj[1,:]
    z2_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_γ0 = arrks_RRK.u_adj[1,:]
    z2_γ0 = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK3 solution",
            arrks_RRK.u[1,:],
            arrks_RRK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RRK3 entropy production",
            ts_RRK.t,arrks_RRK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RRK3 relaxation parameter",
            ts_RRK.t,arrks_RRK.γ,
            xlabel=L"t",
            ylabel=L"\gamma",
            legend=false)
        display(plot!())

        plot(title="RRK3 adjoint solution",
            z1,z2,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK3 adjoint solution (Δt*-const)",
            z1_dt0,z2_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK3 adjoint solution (γ-const)",
            z1_γ0,z2_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK3",ts_RRK.t)
        write(file,"y1_fwd_RRK3",arrks_RRK.u[1,:])
        write(file,"y2_fwd_RRK3",arrks_RRK.u[2,:])
        write(file,"y1_adj_RRK3",z1)
        write(file,"y2_adj_RRK3",z2)
        write(file,"y1_adj_RRK3_dt0",z1_dt0)
        write(file,"y2_adj_RRK3_dt0",z2_dt0)
        write(file,"y1_adj_RRK3_g0",z1_γ0)
        write(file,"y2_adj_RRK3_g0",z2_γ0)
    end


    # RK2 case ----------------------------------------------------------------#
    rk = rk2
    ts_RK.dt = dt_RK2
    ts_RRK.dt = dt_RK2

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)

    if make_plot
        plot(title="RK2 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK2 solution",
            arrks_RK.u[1,:],
            arrks_RK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK2 adjoint solution",
            arrks_RK.u_adj[1,:],
            arrks_RK.u_adj[2,:],
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK2",dt_RK2)
        write(file,"t_RK2",ts_RK.t)
        write(file,"y1_fwd_RK2",arrks_RK.u[1,:])
        write(file,"y2_fwd_RK2",arrks_RK.u[2,:])
        write(file,"y1_adj_RK2",arrks_RK.u_adj[1,:])
        write(file,"y2_adj_RK2",arrks_RK.u_adj[2,:])
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1 = arrks_RRK.u_adj[1,:]
    z2 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_dt0 = arrks_RRK.u_adj[1,:]
    z2_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_γ0 = arrks_RRK.u_adj[1,:]
    z2_γ0 = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK2 solution",
            arrks_RRK.u[1,:],
            arrks_RRK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RRK2 entropy production",
            ts_RRK.t,arrks_RRK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RRK2 relaxation parameter",
            ts_RRK.t,arrks_RRK.γ,
            xlabel=L"t",
            ylabel=L"\gamma",
            legend=false)
        display(plot!())

        plot(title="RRK2 adjoint solution",
            z1,z2,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK2 adjoint solution (Δt*-const)",
            z1_dt0,z2_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK2 adjoint solution (γ-const)",
            z1_γ0,z2_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK2",ts_RRK.t)
        write(file,"y1_fwd_RRK2",arrks_RRK.u[1,:])
        write(file,"y2_fwd_RRK2",arrks_RRK.u[2,:])
        write(file,"y1_adj_RRK2",z1)
        write(file,"y2_adj_RRK2",z2)
        write(file,"y1_adj_RRK2_dt0",z1_dt0)
        write(file,"y2_adj_RRK2_dt0",z2_dt0)
        write(file,"y1_adj_RRK2_g0",z1_γ0)
        write(file,"y2_adj_RRK2_g0",z2_γ0)
    end

    # RK4, big dt -------------------------------------------------------------#
    rk = rk4
    dt_big = 0.9
    ts_RK.dt = dt_big
    ts_RRK.dt = dt_big

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)

    if make_plot
        plot(title="RK4 entropy production (dt big)",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK4 solution (dt big)",
            arrks_RK.u[1,:],
            arrks_RK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK4 adjoint solution (dt big)",
            arrks_RK.u_adj[1,:],
            arrks_RK.u_adj[2,:],
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_big",dt_big)
        write(file,"t_RK_big",ts_RK.t)
        write(file,"y1_fwd_RK_big",arrks_RK.u[1,:])
        write(file,"y2_fwd_RK_big",arrks_RK.u[2,:])
        write(file,"y1_adj_RK_big",arrks_RK.u_adj[1,:])
        write(file,"y2_adj_RK_big",arrks_RK.u_adj[2,:])
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1 = arrks_RRK.u_adj[1,:]
    z2 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_dt0 = arrks_RRK.u_adj[1,:]
    z2_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_γ0 = arrks_RRK.u_adj[1,:]
    z2_γ0 = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK4 solution (dt big)",
            arrks_RRK.u[1,:],
            arrks_RRK.u[2,:],
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 entropy production (dt big)",
            ts_RRK.t,arrks_RRK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RRK4 relaxation parameter (dt big)",
            ts_RRK.t,arrks_RRK.γ,
            xlabel=L"t",
            ylabel=L"\gamma",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (dt big)",
            z1,z2,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (Δt*-const, dt big)",
            z1_dt0,z2_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (γ-const, dt big)",
            z1_γ0,z2_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK_big",ts_RRK.t)
        write(file,"y1_fwd_RRK_big",arrks_RRK.u[1,:])
        write(file,"y2_fwd_RRK_big",arrks_RRK.u[2,:])
        write(file,"y1_adj_RRK_big",z1)
        write(file,"y2_adj_RRK_big",z2)
        write(file,"y1_adj_RRK_dt0_big",z1_dt0)
        write(file,"y2_adj_RRK_dt0_big",z2_dt0)
        write(file,"y1_adj_RRK_g0_big",z1_γ0)
        write(file,"y2_adj_RRK_g0_big",z2_γ0)

        write(file,"eta_RK_big",arrks_RK.Δη)
        write(file,"eta_RRK_big",arrks_RRK.Δη)
        write(file,"gamma_big",arrks_RRK.γ)
        close(file)
    end

    println("End of adjoint growth experiment.")
end
