# This code runs some numerical experiments and generates .mat files for adjoint
# RRK paper. In particular, we run two experiments for a simple nonlinear
# pendulum problem (Ranocha 2019):
#   1. Growth of adjoint solution using RRK w/ & w/o proper linearization.
#   2. Derivative tests for verifying linearization of RRK and quantifying
#      errors from improper linearization.

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

# Flags
write_mat = true  # Want to output .mat files?
run_adj_growth = true # Want to run adjoint growth test?
run_derv_test = true # Want to run derivative tests?
make_plot = false # Want to output plots?

if write_mat
    file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK paper/figs/RRK_nlpen.mat","w")
end

t0 = 0
u0 = [1.5,1]

arrks_RK = AdjRRK_struct()
@pack! arrks_RK = f,df,η,∇η,Hη
@pack! arrks_RK = u0
arrks_RK.return_time = true
arrks_RK.return_Δη = true

ts_RK = Time_struct()
@pack! ts_RK = t0


arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη
@pack! arrks_RRK = u0
arrks_RRK.return_time = true
arrks_RRK.return_Δη = true

ts_RRK = Time_struct()
@pack! ts_RRK = t0



#==============================================================================#
# Running RK and RRK to observe growth of adjoint solution.

ts_RK.T  = 200
ts_RRK.T = 200

dt_RK4 = 0.1#0.9
dt_RK3 = 0.1#0.65
dt_RK2 = 0.1#0.4

if run_adj_growth
    println("Running adjoint growth numerical experiment.")


    # RK4 case ----------------------------------------------------------------#
    rk = rk4
    ts_RK.dt = dt_RK4
    ts_RRK.dt = dt_RK4

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)
    u1_RK = arrks_RK.u[1,:]
    u2_RK = arrks_RK.u[2,:]

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    z1_RK = arrks_RK.u_adj[1,:]
    z2_RK = arrks_RK.u_adj[2,:]

    if make_plot
        plot(title="RK4 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK4 solution",
            u1_RK,u2_RK,
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK4 adjoint solution",
            z1_RK,z2_RK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK4",dt_RK4)
        write(file,"t_RK4",ts_RK.t)
        write(file,"y1_RK4",u1_RK)
        write(file,"y2_RK4",u2_RK)
        write(file,"y1_RK4_adj",z1_RK)
        write(file,"y2_RK4_adj",z2_RK)
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    u1_RRK = arrks_RRK.u[1,:]
    u2_RRK = arrks_RRK.u[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_γ0 = arrks_RRK.u_adj[1,:]
    z2_RRK_γ0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_dt0 = arrks_RRK.u_adj[1,:]
    z2_RRK_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK = arrks_RRK.u_adj[1,:]
    z2_RRK = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK4 solution",
            u1_RRK,u2_RRK,
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

        plot(title="RRK4 adjoint solution (γ-const)",
            z1_RRK_γ0,z2_RRK_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution (Δt*-const)",
            z1_RRK_dt0,z2_RRK_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK4 adjoint solution",
            z1_RRK,z2_RRK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK4",ts_RRK.t)
        write(file,"y1_RRK4",u1_RRK)
        write(file,"y2_RRK4",u2_RRK)
        write(file,"y1_RRK4_adj",z1_RRK)
        write(file,"y2_RRK4_adj",z2_RRK)
        write(file,"y1_RRK4_adj_g0",z1_RRK_γ0)
        write(file,"y2_RRK4_adj_g0",z2_RRK_γ0)
        write(file,"y1_RRK4_adj_dt0",z1_RRK_dt0)
        write(file,"y2_RRK4_adj_dt0",z2_RRK_dt0)
    end


    # RK3 case ----------------------------------------------------------------#
    rk = rk3
    ts_RK.dt = dt_RK3
    ts_RRK.dt = dt_RK3

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)
    u1_RK = arrks_RK.u[1,:]
    u2_RK = arrks_RK.u[2,:]

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    z1_RK = arrks_RK.u_adj[1,:]
    z2_RK = arrks_RK.u_adj[2,:]

    if make_plot
        plot(title="RK3 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK3 solution",
            u1_RK,u2_RK,
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK3 adjoint solution",
            z1_RK,z2_RK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK3",dt_RK3)
        write(file,"t_RK3",ts_RK.t)
        write(file,"y1_RK3",u1_RK)
        write(file,"y2_RK3",u2_RK)
        write(file,"y1_RK3_adj",z1_RK)
        write(file,"y2_RK3_adj",z2_RK)
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    u1_RRK = arrks_RRK.u[1,:]
    u2_RRK = arrks_RRK.u[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_γ0 = arrks_RRK.u_adj[1,:]
    z2_RRK_γ0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_dt0 = arrks_RRK.u_adj[1,:]
    z2_RRK_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK = arrks_RRK.u_adj[1,:]
    z2_RRK = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK3 solution",
            u1_RRK,u2_RRK,
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

        plot(title="RRK3 adjoint solution (γ-const)",
            z1_RRK_γ0,z2_RRK_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK3 adjoint solution (Δt*-const)",
            z1_RRK_dt0,z2_RRK_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK3 adjoint solution",
            z1_RRK,z2_RRK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK3",ts_RRK.t)
        write(file,"y1_RRK3",u1_RRK)
        write(file,"y2_RRK3",u2_RRK)
        write(file,"y1_RRK3_adj",z1_RRK)
        write(file,"y2_RRK3_adj",z2_RRK)
        write(file,"y1_RRK3_adj_g0",z1_RRK_γ0)
        write(file,"y2_RRK3_adj_g0",z2_RRK_γ0)
        write(file,"y1_RRK3_adj_dt0",z1_RRK_dt0)
        write(file,"y2_RRK3_adj_dt0",z2_RRK_dt0)
    end


    # RK2 case ----------------------------------------------------------------#
    rk = rk2
    ts_RK.dt = dt_RK2
    ts_RRK.dt = dt_RK2

    # Running RK
    RK_solver!(arrks_RK,ts_RK,rk)
    u1_RK = arrks_RK.u[1,:]
    u2_RK = arrks_RK.u[2,:]

    # Running adj RK
    arrks_RK.uT_adj = arrks_RK.u[:,end]
    RK_solver!(arrks_RK,ts_RK,rk;adj=true)
    z1_RK = arrks_RK.u_adj[1,:]
    z2_RK = arrks_RK.u_adj[2,:]

    if make_plot
        plot(title="RK2 entropy production",
            ts_RK.t,arrks_RK.Δη,
            xlabel=L"t",
            ylabel=L"\eta(t)-\eta(0)",
            legend=false)
        display(plot!())

        plot(title="RK2 solution",
            u1_RK,u2_RK,
            xlabel=L"y_1",
            ylabel=L"y_2",
            legend=false)
        display(plot!())

        plot(title="RK2 adjoint solution",
            z1_RK,z2_RK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"dt_RK2",dt_RK2)
        write(file,"t_RK2",ts_RK.t)
        write(file,"y1_RK2",u1_RK)
        write(file,"y2_RK2",u2_RK)
        write(file,"y1_RK2_adj",z1_RK)
        write(file,"y2_RK2_adj",z2_RK)
    end

    # Running RRK
    RRK_solver!(arrks_RRK,ts_RRK,rk)
    u1_RRK = arrks_RRK.u[1,:]
    u2_RRK = arrks_RRK.u[2,:]

    # Running adj RRK (γ-constant)
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_γ0 = arrks_RRK.u_adj[1,:]
    z2_RRK_γ0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK (Δt*-constant)
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK_dt0 = arrks_RRK.u_adj[1,:]
    z2_RRK_dt0 = arrks_RRK.u_adj[2,:]

    # Running adj RRK
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    arrks_RRK.uT_adj = arrks_RRK.u[:,end]
    RRK_solver!(arrks_RRK,ts_RRK,rk;adj=true)
    z1_RRK = arrks_RRK.u_adj[1,:]
    z2_RRK = arrks_RRK.u_adj[2,:]

    if make_plot
        plot(title="RRK2 solution",
            u1_RRK,u2_RRK,
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

        plot(title="RRK2 adjoint solution (γ-const)",
            z1_RRK_γ0,z2_RRK_γ0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK2 adjoint solution (Δt*-const)",
            z1_RRK_dt0,z2_RRK_dt0,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())

        plot(title="RRK2 adjoint solution",
            z1_RRK,z2_RRK,
            xlabel=L"\lambda_1",
            ylabel=L"\lambda_2",
            legend=false)
        display(plot!())
    end

    if write_mat
        write(file,"t_RRK2",ts_RRK.t)
        write(file,"y1_RRK2",u1_RRK)
        write(file,"y2_RRK2",u2_RRK)
        write(file,"y1_RRK2_adj",z1_RRK)
        write(file,"y2_RRK2_adj",z2_RRK)
        write(file,"y1_RRK2_adj_g0",z1_RRK_γ0)
        write(file,"y2_RRK2_adj_g0",z2_RRK_γ0)
        write(file,"y1_RRK2_adj_dt0",z1_RRK_dt0)
        write(file,"y2_RRK2_adj_dt0",z2_RRK_dt0)
    end
    println("End of adjoint growth numerical experiment.")
end


#==============================================================================#
# Running derivative test of different lin RRK to observe effects of improper
# linearization.

ts_RK.T  = 200
ts_RRK.T = 200

dt_RK4 = 0.1#0.9
dt_RK3 = 0.1#0.65
dt_RK2 = 0.1#0.4

if run_derv_test
    println("Running derivative numerical tests.")

    Nref = 5
    h0 = 2^(-12)
    if write_mat
        write(file,"Nref",Nref)
        write(file,"h0",h0)
    end

    arrks_h = AdjRRK_struct()
    @pack! arrks_h = f,df,η,∇η,Hη

    Random.seed!(123)
    arrks_RRK.u0_lin = randn(2)


    # RK4 case ----------------------------------------------------------------#
    rk = rk4
    ts_RRK.dt = dt_RK4

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["γ constant" "Δt* constant" "proper lin."]
        markers = [:circle :square :star5]

        errors = [errs_RRK_γ0,errs_RRK_dt0,errs_RRK]
        plot(title="FD error (RRK4)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK_γ0,rate_RRK_dt0,rate_RRK]
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
        write(file,"errs_RRK4_g0",errs_RRK_γ0)
        write(file,"errs_RRK4_dt0",errs_RRK_dt0)
        write(file,"errs_RRK4",errs_RRK)
    end


    # RK3 case ----------------------------------------------------------------#
    rk = rk3
    ts_RRK.dt = dt_RK3

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["γ constant" "Δt* constant" "proper lin."]
        markers = [:circle :square :star5]

        errors = [errs_RRK_γ0,errs_RRK_dt0,errs_RRK]
        plot(title="FD error (RK3)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK_γ0,rate_RRK_dt0,rate_RRK]
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
        write(file,"errs_RRK3_g0",errs_RRK_γ0)
        write(file,"errs_RRK3_dt0",errs_RRK_dt0)
        write(file,"errs_RRK3",errs_RRK)
    end


    # RK2 case ----------------------------------------------------------------#
    rk = rk2
    ts_RRK.dt = dt_RK2

    # γ-const case
    arrks_RRK.γ_cnst = true
    arrks_RRK.dt_cnst = true
    errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # dt*-const case
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = true
    errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    # with proper linearization
    arrks_RRK.γ_cnst = false
    arrks_RRK.dt_cnst = false
    errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk,h0,Nref)

    if make_plot
        labels = ["γ constant" "Δt* constant" "proper lin."]
        markers = [:circle :square :star5]

        errors = [errs_RRK_γ0,errs_RRK_dt0,errs_RRK]
        plot(title="FD error (RK2)",
            h,errors,
            xlabel=L"h",
            ylabel="error",
            label=labels,
            marker=markers,
            xaxis=:log,
            yaxis=:log)
        display(plot!())

        # rates = [rate_RRK_γ0,rate_RRK_dt0,rate_RRK]
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
        write(file,"errs_RRK2_g0",errs_RRK_γ0)
        write(file,"errs_RRK2_dt0",errs_RRK_dt0)
        write(file,"errs_RRK2",errs_RRK)
    end

    if write_mat
        close(file)
    end
    println("End of derivative numerical tests.")
end
