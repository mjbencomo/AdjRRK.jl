## Algorithms for simple explicit RK methods
# function RK_update(u,f,dt,rk::RKs)
#     U = u
#     F = f(U)
#     du = rk.b[1].*F
#
#     for s=2:rk.stages
#         U = @. u + dt*rk.a[s-1]*F
#         F = f(U)
#         du = @. du + rk.b[s]*F
#     end
#     return @. u + dt*du
# end
#
# function RK_update_lin(flds,ops,dt,rk::RKs)
#     w,u  = flds
#     f,df = ops
#     W = w
#     U = u
#     JW = df(U,W)
#     dw = rk.b[1].*JW
#     for s=2:rk.stages
#         W = @. w + dt*rk.a[s-1]*JW
#         F = f(U)
#         U = @. u + dt*rk.a[s-1]*F
#         JW = df(U,W)
#         dw = @. dw + rk.b[s]*JW
#     end
#     return @. w + dt*dw
# end
#
# function RK_update_adj(flds,ops,dt,rk::RKs)
#     z,u  = flds
#     f,df = ops
#     a = rk.a
#     b = rk.b
#     S = rk.stages
#
#     # recomputing fwd internal stages
#     U = zeros(length(u),S)
#     U[:,1] = u
#     for s=2:S
#         U[:,s] =  u .+ dt.*a[s-1].*f(U[:,s-1])
#     end
#
#     Z = df(U[:,S],z;adj=true)
#     Z = @. dt*b[S]*Z
#     dz = Z
#     for s=S-1:-1:1
#         Z = @. b[s]*z + a[s]*Z
#         Z = df(U[:,s],Z;adj=true)
#         Z = @. dt*Z
#         dz = @. dz + Z
#     end
#     return @. z + dz
# end


## Algorithms of generic explicit RK methods
function RK_update(u,f,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    U = zeros(length(u),s)
    U[:,1] = u
    du = b[1].*f(U[:,1])

    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f(U[:,j])
        end
        U[:,i] = u + dt*U[:,i]
        du += b[i]*f(U[:,i])
    end
    return @. u + dt*du
end

function RK_update_lin(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u  = flds
    f,df = ops

    #recomputing internal fwd stages
    U = zeros(length(u),s)
    U[:,1] = u
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f(U[:,j])
        end
        U[:,i] = u + dt*U[:,i]
    end

    #linear internal stages
    W = zeros(length(u),s)
    W[:,1] = w
    dw = b[1].*df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*df(U[:,j],W[:,j])
        end
        W[:,i] = w + dt*W[:,i]
        dw += b[i]*df(U[:,i],W[:,i])
    end
    return @. w + dt*dw
end

function RK_update_adj(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    z,u  = flds
    f,df = ops

    #recomputing fwd internal stages
    U = zeros(length(u),s)
    U[:,1] = u
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f(U[:,j])
        end
        U[:,i] = u + dt*U[:,i]
    end

    #adjoint internal stages
    Z = zeros(length(u),s)
    Z[:,s] = dt*b[s]*df(U[:,s],z;adj=true)
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = b[i]*z + Z[:,i]
        Z[:,i] = dt*df(U[:,i],Z[:,i];adj=true)
        dz += Z[:,i]
    end
    return @. z + dz
end


## RK SOLVER
#
# INPUTS:
#   arrks::AdjRRK_struct
#       [RK fwd case]
#           u0,     initial condition
#           f,      RHS function
#           η,      entropy function (if return_Δη=true)
#           return_time,
#           return_Δη,
#       [RK lin case; if lin=true]
#           u,      solution of forward RK
#           u0_lin, initial condition for linearized solution
#           f,      RHS function
#           df,     Jacobian of RHS function
#       [RK adj case; if adj=true]
#           u,      solution of forward RK
#           uT_adj, final condition for adjoint solution
#           f,      RHS function
#           df,     Jacobian of RHS function
#   ts::Time_struct
#       [RK fwd case]
#           t0,     intial time
#           T,      final time
#           dt,     time step size
#       [RK lin case; if lin=true]
#           dt,     time step size
#           Nt,     number of time steps
#       [RK adj case; if adj=true]
#           dt,     time step size
#           Nt,     number of time steps
#   rk::RKs
#
# OUTPUTS:
#   arrks::AdjRRK_struct
#       [RK fwd case]
#           u,      solution
#           Δη,     change in entropy (if return_Δη=true)
#       [RK lin case; if lin=true]
#           u_lin,  linearized solution
#       [RK adj case; if adj=true]
#           u_adj,  adjoint solution
#   ts::Time_struct
#       [RK fwd case]
#           dt,     may get modified
#           Nt,     number of time steps
#           t,      time grid points (if return_time=true)

function RK_solver!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct;lin=false,adj=false)
    if lin && ~(adj)
        RK_lin!(arrks,ts,rk)
    elseif adj
        RK_adj!(arrks,ts,rk)
    else
        RK_fwd!(arrks,ts,rk)
    end
end

function RK_fwd!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack t0,T,dt = ts
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)
    @pack! ts = dt,Nt

    @unpack f = arrks
    @unpack u0 = arrks

    u = zeros(length(u0),Nt)
    u[:,1] = u0

    for k=1:Nt-1
        u[:,k+1] = RK_update(u[:,k],f,dt,rk)
    end
    @pack! arrks = u

    if arrks.return_time
        t = range(t0,stop=T,length=Nt) |> collect
        @pack! ts = t
    end
    if arrks.return_Δη
        @unpack η = arrks
        Δη = zeros(Nt)
        η0 = η(u[:,1])
        for k=1:Nt-1
            Δη[k+1] = η(u[:,k+1])-η0
        end
        @pack! arrks = Δη
    end
end

function RK_lin!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack u0_lin,u = arrks

    u_lin = zeros(size(u))
    u_lin[:,1] = u0_lin

    for k=1:Nt-1
        flds = (u_lin[:,k],u[:,k])
        u_lin[:,k+1] = RK_update_lin(flds,(f,df),dt,rk)
    end
    @pack! arrks = u_lin
end

function RK_adj!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack uT_adj,u = arrks

    u_adj = zeros(size(u))
    u_adj[:,Nt] = uT_adj

    for k=Nt-1:-1:1
        flds = (u_adj[:,k+1],u[:,k])
        u_adj[:,k] = RK_update_adj(flds,(f,df),dt,rk)
    end
    @pack! arrks = u_adj
end
