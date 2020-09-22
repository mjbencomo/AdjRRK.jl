# Algorithms for simple explicit RK methods.

function RK_update(u,f,dt,rk::RKs)
    U = u
    F = f(U)
    du = rk.b[1].*F

    for s=2:rk.stages
        U = @. u + dt*rk.a[s-1]*F
        F = f(U)
        du = @. du + rk.b[s]*F
    end
    return @. u + dt*du
end

function RK_update_lin(flds,ops,dt,rk::RKs)
    w,u  = flds
    f,df = ops
    W = w
    U = u
    JW = df(U,W)
    dw = rk.b[1].*JW
    for s=2:rk.stages
        W = @. w + dt*rk.a[s-1]*JW
        F = f(U)
        U = @. u + dt*rk.a[s-1]*F
        JW = df(U,W)
        dw = @. dw + rk.b[s]*JW
    end
    return @. w + dt*dw
end

function RK_update_adj(flds,ops,dt,rk::RKs)
    z,u  = flds
    f,df = ops
    a = rk.a
    b = rk.b
    S = rk.stages

    # recomputing fwd internal stages
    U = zeros(length(u),S)
    U[:,1] = u
    for s=2:S
        U[:,s] =  u .+ dt.*a[s-1].*f(U[:,s-1])
    end

    Z = df(U[:,S],z;adj=true)
    Z = @. dt*b[S]*Z
    dz = Z
    for s=S-1:-1:1
        Z = @. b[s]*z + a[s]*Z
        Z = df(U[:,s],Z;adj=true)
        Z = @. dt*Z
        dz = @. dz + Z
    end
    return @. z + dz
end



# --------------------------------------#
# RK ALGORITHM
#
# Remark: in_flds and ops depend on the optional flags.
#
# INPUT ARGS:
#   in_flds = u0,
#           = (w0,u), [lin=true]
#           = (zT,u), [adj=true]
#       ops = f,
#           = (f,η),  [return_Δη=true]
#           = (f,df), [lin=true or adj=true]
#      Time = (t0,T,dt), time axis info.
#       rk = RKs struct storing coefficients.
# OPTIONAL FLAGS:
#         lin = running linearized algorithm
#         adj = running adjoint algorithm
# return_time = return time axis
#   return_Δη = return change in entropy
#
# OUTPUT:
#   u,
#   (u,Δu), return_Δη=true
#   (u,t),  return_time=true
#   w,      lin=true
#   z,      adj=true
# ---------------------------------------#
function RK_solver(in_flds,ops,Time,rk;
lin=false,
adj=false,
return_time=false,
return_Δη=false)

    t0,T,dt = Time
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)

    if lin && ~(adj)
    # running linearized RK
        f,df = ops
        w0,u = in_flds
        w = zeros(size(u))
        w[:,1] = w0
        for k=1:Nt-1
            flds = (w[:,k],u[:,k])
            w[:,k+1] = RK_update_lin(flds,(f,df),dt,rk)
        end
        out_flds = w

    elseif adj
    # running adjoint RK
        f,df = ops
        zT,u = in_flds
        z = zeros(size(u))
        z[:,Nt] = zT
        for k=Nt-1:-1:1
            flds = (z[:,k+1],u[:,k])
            z[:,k] = RK_update_adj(flds,(f,df),dt,rk)
        end
        out_flds = z

    else
    # running standard RK
        u0 = in_flds
        u = zeros(length(u0),Nt)
        u[:,1] = u0
        if return_Δη #case where Δη is tracked
            f,η = ops
            Δηu = zeros(Nt)
            ηu0 = η(u[:,1])
            for k=1:Nt-1
                u[:,k+1] = RK_update(u[:,k],f,dt,rk)
                Δηu[k+1] = η(u[:,k+1])-ηu0
            end
            out_flds = (u,Δηu)

            if return_time
                t = range(t0,stop=T,length=Nt) |> collect
                out_flds = (out_flds...,t)
            end

        else #case where Δη is not tracked
            f = ops
            for k=1:Nt-1
                u[:,k+1] = RK_update(u[:,k],f,dt,rk)
            end
            out_flds = u

            if return_time
                t = range(t0,stop=T,length=Nt) |> collect
                out_flds = (out_flds,t)
            end
        end
    end
    return out_flds
end
