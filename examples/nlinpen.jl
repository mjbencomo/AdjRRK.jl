# Testing out RRK method (linearized and adjoint code) on a small nonlinear ODE
# system (nlinear pendulum) from Ranocha 2019.

# Comparing the linearized and adjoint solutions from RRK algs when relaxation
# parameter γ is viewed as a constant or function of state fields during
# linearization.

using AdjRRK
using Plots

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


rk = rk4
u0 = [1.5,1]

t0 = 0
T  = 500


## Computing high res solution with standard RK.
dt = 0.005
Time = (t0,T,dt)

u,Δη,t = RK_solver(u0,(f,η),Time,rk; return_Δη=true,return_time=true)
w      = RK_solver((u0,u),(f,df),Time,rk; lin=true)
z      = RK_solver((u[:,end],u),(f,df),Time,rk; adj=true)

# display(plot(t,Δη,
#     title="Change in entropy",
#     label="RK (dt=$dt)",
#     xlabel="t",
#     ylabel="Δη"))
#
# display(plot(u[1,:],u[2,:],
#     title="solution",
#     label="RK (dt=$dt)",
#     xlabel="u1",
#     ylabel="u2"))
#
# display(plot(w[1,:],w[2,:],
#     title="linearized solution",
#     label="RK (dt=$dt)",
#     xlabel="w1",
#     ylabel="w2"))
#
# display(plot(z[1,:],z[2,:],
#     title="adjoint solution",
#     label="RK (dt=$dt)",
#     xlabel="z1",
#     ylabel="z2"))


## Computing low res solution with standard RK.
dt_RK = 0.5
Time = (t0,T,dt_RK)

u_RK,Δη_RK,t_RK = RK_solver(u0,(f,η),Time,rk; return_Δη=true,return_time=true)
w_RK            = RK_solver((u0,u_RK),(f,df),Time,rk; lin=true)
z_RK            = RK_solver((u_RK[:,end],u_RK),(f,df),Time,rk; adj=true)


## Computing solution using IDT.

u_IDT,γ_IDT,Δη_IDT = IDT_solver(u0,(f,η,∇η),Time,rk; return_Δη=true)
w_IDTγ0 = IDT_solver((u0,u_IDT,γ_IDT),(f,df),Time,rk; lin=true,γcnst=true)
w_IDT   = IDT_solver((u0,u_IDT,γ_IDT),(f,df,∇η,Hη),Time,rk; lin=true)
z_IDTγ0 = IDT_solver((u_IDT[:,end],u_IDT,γ_IDT),(f,df),Time,rk; adj=true,γcnst=true)
z_IDT   = IDT_solver((u_IDT[:,end],u_IDT,γ_IDT),(f,df,∇η,Hη),Time,rk; adj=true)

## Computing solution using RRK.

u_RRK,γ_RRK,Δη_RRK,t_RRK,dt_corr = RRK_solver(u0,(f,η,∇η),Time,rk; return_Δη=true, return_time=true)
w_RRKγ0 = RRK_solver((u0,u_RRK,γ_RRK,dt_corr),(f,df),Time,rk; lin=true,γcnst=true)
w_RRK   = RRK_solver((u0,u_RRK,γ_RRK,dt_corr),(f,df,∇η,Hη),Time,rk; lin=true)
z_RRKγ0 = RRK_solver((u_RRK[:,end],u_RRK,γ_RRK,dt_corr),(f,df),Time,rk; adj=true,γcnst=true)
z_RRK   = RRK_solver((u_RRK[:,end],u_RRK,γ_RRK,dt_corr),(f,df,∇η,Hη),Time,rk; adj=true)


## PLOTS

#plotting change in entropy
plot(title="Change in entropy",
    xlabel="t",
    ylabel="Δη")

plot!(t,Δη,label="RK (dt=$dt)")
# plot!( t_RK,  Δη_RK, label="RK (dt=$dt_RK)")
plot!( t_RK, Δη_IDT, label="IDT")
plot!(t_RRK, Δη_RRK, label="RRK")
display(plot!())
#png("entropy")

#plotting solutions
plot(title="solution",
    xlabel="u1",
    ylabel="u2")
plot!(    u[1,:],    u[2,:], label="RK (dt=$dt)")
plot!( u_RK[1,:], u_RK[2,:], label="RK (dt=$dt_RK)")
display(plot!())
# png("u_RK")

plot(title="solution",
    xlabel="u1",
    ylabel="u2")
plot!(    u[1,:],    u[2,:], label="RK (dt=$dt)")
plot!(u_IDT[1,:],u_IDT[2,:], label="IDT")
display(plot!())
# png("u_IDT")

plot(title="solution",
    xlabel="u1",
    ylabel="u2")
plot!(    u[1,:],    u[2,:], label="RK (dt=$dt)")
plot!(u_RRK[1,:],u_RRK[2,:], label="RRK")
display(plot!())
# png("u_RRK")


#plotting linearized solution
plot(title="linearized solution",
    xlabel="w1",
    ylabel="w2")
plot!(      w[1,:],      w[2,:], label="RK (dt=$dt)")
plot!(   w_RK[1,:],   w_RK[2,:], label="RK (dt=$dt_RK)")
display(plot!())
# png("w_RK")

plot(title="linearized solution",
    xlabel="w1",
    ylabel="w2")
plot!(      w[1,:],      w[2,:], label="RK (dt=$dt)")
plot!(w_IDTγ0[1,:],w_IDTγ0[2,:], label="IDTγ0")
display(plot!())
# png("w_IDTcnst")

plot(title="linearized solution",
    xlabel="w1",
    ylabel="w2")
plot!(      w[1,:],      w[2,:], label="RK (dt=$dt)")
plot!(w_RRKγ0[1,:],w_RRKγ0[2,:], label="RRKγ0")
display(plot!())
# png("w_RRKcnst")

plot(title="linearized solution",
    xlabel="w1",
    ylabel="w2")
plot!(      w[1,:],      w[2,:], label="RK (dt=$dt)")
plot!(  w_IDT[1,:],  w_IDT[2,:], label="IDT")
display(plot!())
# png("w_IDT")

plot(title="linearized solution",
    xlabel="w1",
    ylabel="w2")
plot!(      w[1,:],      w[2,:], label="RK (dt=$dt)")
plot!(  w_RRK[1,:],  w_RRK[2,:], label="RRK")
display(plot!())
# png("w_RRK")


#plotting adjoint solution
plot(title="adjoint solution",
    xlabel="z1",
    ylabel="z2")
plot!(      z[1,:],      z[2,:], label="RK (dt=$dt)")
plot!(   z_RK[1,:],   z_RK[2,:], label="RK (dt=$dt_RK)")
display(plot!())
# png("z_RK")

plot(title="adjoint solution",
    xlabel="z1",
    ylabel="z2")
plot!(      z[1,:],      z[2,:], label="RK (dt=$dt)")
plot!(z_IDTγ0[1,:],z_IDTγ0[2,:], label="IDTγ0")
display(plot!())
# png("z_IDTcnst")

plot(title="adjoint solution",
    xlabel="z1",
    ylabel="z2")
plot!(      z[1,:],      z[2,:], label="RK (dt=$dt)")
plot!(z_RRKγ0[1,:],z_RRKγ0[2,:], label="RRKγ0")
display(plot!())
# png("z_RRKcnst")

plot(title="adjoint solution",
    xlabel="z1",
    ylabel="z2")
plot!(      z[1,:],      z[2,:], label="RK (dt=$dt)")
plot!(  z_IDT[1,:],  z_IDT[2,:], label="IDT")
display(plot!())
# png("z_IDT")

plot(title="adjoint solution",
    xlabel="z1",
    ylabel="z2")
plot!(      z[1,:],      z[2,:], label="RK (dt=$dt)")
plot!(  z_RRK[1,:],  z_RRK[2,:], label="RRK")
display(plot!())
#png("z_RRK")


# display(plot(t,γ,title="relaxation parameter",legend=false))
#
# plot(u_true[1,:],u_true[2,:],title="nonlinear solution",label="RK (dt=0.05)")
# plot!(u_RK[1,:],u_RK[2,:],label="RK")
# display(plot!(u_RRK[1,:],u_RRK[2,:],label="RRK"))
#
# plot(w_true[1,:],w_true[2,:],title="linearized solution",label="RK (dt=0.05)")
# plot!(w_RK[1,:],w_RK[2,:],label="RK")
# plot!(w_RRKγ0[1,:],w_RRKγ0[2,:],label="RRKγ0")
# display(plot!(w_RRK[1,:],w_RRK[2,:],label="RRK"))
#
# plot(z_true[1,:],z_true[2,:],title="adjoint solution",label="RK (dt=0.05)")
# plot!(z_RK[1,:],z_RK[2,:],label="RK")
# plot!(z_RRKγ0[1,:],z_RRKγ0[2,:],label="RRKγ0")
# display(plot!(z_RRK[1,:],z_RRK[2,:],label="RRK"))
#
#
# plot(w_true[1,:],w_true[2,:],title="comparing linearized solutions",label="RK (dt=0.05)")
# plot!(w_RK[1,:],w_RK[2,:],label="RK")
# display(plot!(w_RRK[1,:],w_RRK[2,:],label="RRK"))
#
# plot(z_true[1,:],z_true[2,:],title="comparing adjoint solutions",label="RK (dt=0.05)")
# plot!(z_RK[1,:],z_RK[2,:],label="RK")
# display(plot!(z_RRK[1,:],z_RRK[2,:],label="RRK"))
