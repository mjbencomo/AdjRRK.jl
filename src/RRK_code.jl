# Algorithms for simple explicit RRK methods.

function RK_int(u,f,dt,rk::RKs)
    S = rk.stages
    U = zeros(length(u),S)
    F = zeros(length(u),S)

    U[:,1] = u
    F[:,1] = f(u)
    for s=2:S
        U[:,s] = @. u + dt*rk.a[s-1]*F[:,s-1]
        F[:,s] = f(U[:,s])
    end
    return U,F
end

function RK_int_lin(flds,ops,dt,rks::RKs)
    w,u  = flds
    f,df = ops
    S = rks.stages
    U = zeros(length(u),S)
    F = zeros(length(u),S)
    W = zeros(length(u),S)
    JW = zeros(length(u),S)

    U[:,1] = u
    F[:,1] = f(U[:,1])
    W[:,1] = w
    JW[:,1] = df(U[:,1],W[:,1])
    for s=2:S
        U[:,s] = @. u + dt*rks.a[s-1]*F[:,s-1]
        F[:,s] = f(U[:,s])

        W[:,s] = @. w + dt*rks.a[s-1]*JW[:,s-1]
        JW[:,s] = df(U[:,s],W[:,s])
    end
    return U,F,W,JW
end

function bisect(γL,γR,r;rtol=1e-15,γtol=1e-15,nmax=10000)
    for n=1:nmax
        γm = (γL+γR)*.5
        rm = r(γm)

        if abs(rm)<rtol #|| abs((γR-γL)*.5)<γtol
            return γm
        end
        if sign(rm)==sign(r(γL))
            γL = γm
        else
            γR = γm
        end
     end
     error("Exiting from bisect, could not find root under $nmax iterations.\\ γ=$γm, r(γ)=$rm")
end

function RRK_update(flds,ops,dt,rks::RKs)
    u,γ0   = flds
    f,η,∇η = ops
    U,F = RK_int(u,f,dt,rks)

    d = zeros(length(u))
    e = 0
    for s=1:rks.stages
        d = @. d + rk.b[s].*F[:,s]
        e = e + rk.b[s]*( ∇η(U[:,s])⋅F[:,s] )
    end
    d = @. dt*d
    e = dt*e
    ηu = η(u)

    function r(γ)
        return η(@. u + γ*d) - ηu - γ*e
    end

    γ = 1
    #this tol should be relative to η(u0)
    if abs(r(1)) > 1e-15
        γL = 0.9*γ0
        γR = 1.1*γ0
        i = 1
        imax = 100
        while r(γL)>0 && i<=imax
            γR = γL
            γL = γL*0.5
            i  = i+1
            if(i==imax)
                @show γ0
                @show γL
                @show r(γ0)
                @show r(γL)
                @show r(0)
                @show r(0.5)
                @show r(1)
                error("Exiting from RRK_update, could not find γL<0 under 100 iterations")
            end
        end

        i = 1
        while r(γR)<0 && i<=imax
            γL = γR
            γR = γR*2
            i  = i+1
            if i==imax
                error("Exiting from RRK_update, could not find γR>0 under 100 iterations")
            end
        end
        γ = bisect(γL,γR,r)
    end
    u = u + γ*d

    return u,γ
end

function RRKγ0_update_lin(flds,ops,dt,rk::RKs)
    w,u,γ = flds
    new_b = rk.b .* γ
    new_rk = RKs(rk.stages,new_b,rk.a)
    return RK_update_lin((w,u),ops,dt,new_rk)
end

function RRKγ0_update_adj(flds,ops,dt,rk::RKs)
    z,u,γ = flds
    new_b = rk.b .* γ
    new_rk = RKs(rk.stages,new_b,rk.a)
    return RK_update_adj((z,u),ops,dt,new_rk)
end

function RRK_update_lin(flds,ops,dt,rk::RKs)
    w,u,γ,ϱ = flds
    f,df,∇η,Hη = ops
    S = rk.stages
    b = rk.b
    a = rk.a

    u1 = u[:,1]
    u2 = u[:,2]

    U,F,W,JW = RK_int_lin((w,u1),(f,df),dt,rk)

    ∇ηu2 = ∇η(u2)
    Γw = ( ∇ηu2 - ∇η(u1) )⋅w
    Γ̃W = 0
    rγ = 0
    for s=1:S
        ∇ηd = ∇ηu2-∇η(U[:,s])
        rγ = rγ + b[s]*( ∇ηd⋅F[:,s] )
        Γ̃W = Γ̃W + b[s]*( ∇ηd⋅JW[:,s] - F[:,s]⋅Hη(U[:,s],W[:,s]) )
    end
    rγ = dt*rγ
    Γ̃W = γ*dt*Γ̃W
    ρ = Γ̃W + Γw
    ρ = -ρ/rγ

    dw = zeros(length(w))
    for s=1:S
        dw = @. dw + b[s]*( γ*JW[:,s] + ρ*F[:,s] )
    end
    ϱ += ρ
    return (@. w + dt*dw, ϱ )
end

function RRK_update_adj(flds,ops,dt,rk::RKs)
    z,u,γ,ζ = flds
    f,df,∇η,Hη = ops
    S = rk.stages
    b = rk.b
    a = rk.a

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = RK_int(u1,f,dt,rk)

    ∇ηu2 = ∇η(u2)
    Γ̃  = zeros(length(u1),S)
    rγ = 0
    ξ  = 0
    for s=1:S
        ∇ηd    = ∇ηu2-∇η(U[:,s])
        Γ̃[:,s] = b[s].*( df(U[:,s],∇ηd;adj=true) - Hη(U[:,s],F[:,s];adj=true) )
        rγ     = rγ + b[s]*( ∇ηd⋅F[:,s] )
        ξ      = ξ  + b[s]*( F[:,s]⋅z )
    end
    rγ = dt*rγ
    Γ̃  = -γ*dt/rγ .* Γ̃
    Γ  = ( ∇η(u1)-∇ηu2 )./rγ
    ξ  = dt*ξ - ζ

    Z = df(U[:,S],z;adj=true)
    Z = @. dt*γ*b[S]*Z + ξ*Γ̃[:,S]
    dz = Z
    for s=S-1:-1:1
        Z = @. γ*b[s]*z + a[s]*Z
        Z = df(U[:,s],Z;adj=true)
        Z = @. dt*Z + ξ*Γ̃[:,s]
        dz = @. dz + Z
    end
    return @. z + dz + ξ*Γ
end

function RRK_update_lin_last(flds,ops,dt_corr,rk::RKs)
    w,u,γ,ϱ,dt_old = flds
    f,df,∇η,Hη = ops
    S = rk.stages
    b = rk.b
    a = rk.a

    u1 = u[:,1]
    u2 = u[:,2]

    U = zeros(length(u1),S)
    F = zeros(length(u1),S)
    W = zeros(length(u1),S)
    JW = zeros(length(u1),S)

    U[:,1] = u1
    F[:,1] = f(U[:,1])
    W[:,1] = w
    JW[:,1] = df(U[:,1],W[:,1])
    for s=2:S
        U[:,s] = @. u1 + dt_corr*a[s-1]*F[:,s-1]
        F[:,s] = f(U[:,s])
        W[:,s] = @. w + dt_corr*a[s-1]*( JW[:,s-1] - dt_old/dt_corr*ϱ*F[:,s-1] )
        JW[:,s] = df(U[:,s],W[:,s])
    end

    ∇ηu2 = ∇η(u2)
    Γw = ( ∇ηu2 - ∇η(u1) )⋅w
    Γ̃W = 0
    rγ = 0
    for s=1:S
        ∇ηd = ∇ηu2-∇η(U[:,s])
        rγ = rγ + b[s]*( ∇ηd⋅F[:,s] )
        Γ̃W = Γ̃W + b[s]*( ∇ηd⋅JW[:,s] - F[:,s]⋅Hη(U[:,s],W[:,s]) )
    end
    rγ = dt_corr*rγ
    Γ̃W = γ*dt_corr*Γ̃W
    ρ = Γ̃W + Γw
    ρ = -ρ/rγ

    dw = zeros(length(w))
    for s=1:S
        dw = @. dw + b[s]*( γ*JW[:,s] + ρ*F[:,s] )
    end
    return @. w + dt_corr*dw
end

function RRK_update_adj_last(flds,ops,dt_corr,rk::RKs)
    z,u,γ,dt_old = flds
    f,df,∇η,Hη = ops
    S = rk.stages
    b = rk.b
    a = rk.a

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = RK_int(u1,f,dt_corr,rk)

    ∇ηu2 = ∇η(u2)
    Γ̃  = zeros(length(u1),S)
    rγ = 0
    ξ  = 0
    for s=1:S
        ∇ηd    = ∇ηu2-∇η(U[:,s])
        Γ̃[:,s] = b[s].*( df(U[:,s],∇ηd;adj=true) - Hη(U[:,s],F[:,s];adj=true) )
        rγ     = rγ + b[s]*( ∇ηd⋅F[:,s] )
        ξ      = ξ  + b[s]*( F[:,s]⋅z )
    end
    rγ = dt_corr*rγ
    Γ̃  = -γ*dt_corr/rγ .* Γ̃
    Γ  = ( ∇η(u1)-∇ηu2 )./rγ
    ξ  = dt_corr*ξ

    Z = df(U[:,S],z;adj=true)
    Z = @. dt_corr*γ*b[S]*Z + ξ*Γ̃[:,S]
    dz = Z
    ζ = 0
    for s=S-1:-1:1
        ζ += a[s]*( F[:,s]⋅Z )
        Z = @. γ*b[s]*z + a[s]*Z
        Z = df(U[:,s],Z;adj=true)
        Z = @. dt_corr*Z + ξ*Γ̃[:,s]
        dz = @. dz + Z
    end
    ζ *= dt_old
    return (@. z + dz + ξ*Γ, ζ)
end

# --------------------------------------#
# IDT ALGORITHM
#
# Remark: in_flds and ops depend on the optional flags.
#
# INPUT ARGS:
#   in_flds = u0,
#           = (w0,u,γ),     [lin=true]
#           = (zT,u,γ),     [adj=true]
#       ops = (f,η,∇η)
#           = (f,df),       [lin=true or adj=true, & γcnst=true]
#           = (f,df,∇η,Hη), [lin=true or adj=true]
#      Time = (t0,T,dt), time axis info.
#       rk = RKs struct storing coefficients.
# OPTIONAL FLAGS:
#         lin = running linearized algorithm
#         adj = running adjoint algorithm
# return_time = return time axis
#   return_Δη = return change in entropy
#       γcnst = view γ constant during linearization/adjoint
#
# OUTPUT:
#   u,γ
#   (u,γ,Δη), return_Δη=true
#   (u,γ,t),  return_time=true
#   (u,γ,Δη,t)
#   w,      lin=true
#   z,      adj=true
# ---------------------------------------#
function IDT_solver(in_flds,ops,Time,rk;
lin=false,
adj=false,
γcnst=false,
return_time=false,
return_Δη=false)

    t0,T,dt = Time
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)

    if lin && ~(adj)
    #Linearized IDT algorithm
        w0,u,γ = in_flds
        w = zeros(size(u))
        w[:,1] = w0

        if γcnst
            for k=1:Nt-1
                flds = (w[:,k],u[:,k],γ[k+1])
                w[:,k+1] = RRKγ0_update_lin(flds,ops,dt,rk)
            end
            out_flds = w
        else
            ϱ = 0
            for k=1:Nt-1
                flds = (w[:,k],u[:,k:k+1],γ[k+1],ϱ)
                w[:,k+1],ϱ = RRK_update_lin(flds,ops,dt,rk)
            end
            out_flds = w
        end

    elseif adj
    #Adjoint IDT algorithm
        zT,u,γ = in_flds
        z = zeros(size(u))
        z[:,Nt] = zT

        if γcnst
            for k=Nt-1:-1:1
                flds = (z[:,k+1],u[:,k],γ[k+1])
                z[:,k] = RRKγ0_update_adj(flds,ops,dt,rk)
            end
            out_flds = z
        else
            for k=Nt-1:-1:1
                flds = (z[:,k+1],u[:,k:k+1],γ[k+1],0)
                z[:,k] = RRK_update_adj(flds,ops,dt,rk)
            end
            out_flds = z
        end

    else
    #IDT algorithm
        f,η,∇η = ops
        u0     = in_flds
        u      = zeros(length(u0),Nt)
        u[:,1] = u0
        γ      = ones(Nt)
        for k=1:Nt-1
            flds = (u[:,k],γ[k])
            u[:,k+1],γ[k+1] = RRK_update(flds,ops,dt,rk)
        end
        out_flds = (u,γ)

        if return_Δη
            Δηu = zeros(Nt)
            ηu0 = η(u[:,1])
            for k=1:Nt-1
                Δηu[k+1] = η(u[:,k+1])-ηu0
            end
            out_flds = (out_flds...,Δηu)
        end
        if return_time
            t = range(t0,stop=T,length=Nt) |> collect
            out_flds = (out_flds...,t)
        end
    end
    return out_flds
end

# --------------------------------------#
# RRK ALGORITHM
#
# Remark: in_flds and ops depend on the optional flags.
#
# INPUT ARGS:
#   in_flds = u0,
#           = (w0,u,γ,last_dt),     [lin=true]
#           = (zT,u,γ,last_dt),     [adj=true]
#       ops = (f,η,∇η)
#           = (f,df),       [lin=true or adj=true, & γcnst=true]
#           = (f,df,∇η,Hη), [lin=true or adj=true]
#      Time = (t0,T,dt), time axis info.
#       rk = RKs struct storing coefficients.
# OPTIONAL FLAGS:
#         lin = running linearized algorithm
#         adj = running adjoint algorithm
# return_time = return time axis
#   return_Δη = return change in entropy
#       γcnst = view γ constant during linearization/adjoint
#
# OUTPUT:
#   u,γ
#   (u,γ,Δη), return_Δη=true
#   (u,γ,t,last_dt),  return_time=true
#   (u,γ,Δη,t,last_dt)
#   w,      lin=true
#   z,      adj=true
# ---------------------------------------#
function RRK_solver(in_flds,ops,Time,rk;
lin=false,
adj=false,
γcnst=false,
corr=true,
return_time=false,
return_Δη=false)

    t0,T,dt = Time
    Nt_tmp = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt_tmp-1)

    if lin && ~(adj)
    #linearized RRK algorithm
        w0,u,γ,dt_corr = in_flds
        Nt = length(γ)
        w = zeros(size(u))
        w[:,1] = w0

        if γcnst
            for k=1:Nt-2
                flds = (w[:,k],u[:,k],γ[k+1])
                w[:,k+1] = RRKγ0_update_lin(flds,ops,dt,rk)
            end
            #last step
            flds = (w[:,Nt-1],u[:,Nt-1],γ[Nt])
            w[:,Nt] = RRKγ0_update_lin(flds,ops,dt_corr,rk)
            out_flds = w
        else
            ϱ = 0
            for k=1:Nt-2
                flds = (w[:,k],u[:,k:k+1],γ[k+1],ϱ)
                w[:,k+1],ϱ = RRK_update_lin(flds,ops,dt,rk)
            end

            if ~corr
                ϱ = 0
            end
            flds = (w[:,Nt-1],u[:,Nt-1:Nt],γ[Nt],ϱ,dt)
            w[:,Nt] = RRK_update_lin_last(flds,ops,dt_corr,rk)
            out_flds = w
        end

    elseif adj
    #adjoint RRK algorithm
        zT,u,γ,dt_corr = in_flds
        Nt = length(γ)
        z = zeros(size(u))
        z[:,Nt] = zT

        if γcnst
            flds = (z[:,Nt],u[:,Nt-1],γ[Nt])
            z[:,Nt-1] = RRKγ0_update_adj(flds,ops,dt_corr,rk)
            for k=Nt-2:-1:1
                flds = (z[:,k+1],u[:,k],γ[k+1])
                z[:,k] = RRKγ0_update_adj(flds,ops,dt,rk)
            end
            out_flds = z
        else
            flds = (z[:,Nt],u[:,Nt-1:Nt],γ[Nt],dt)
            z[:,Nt-1],ζ = RRK_update_adj_last(flds,ops,dt_corr,rk)

            if ~corr
                ζ = 0
            end
            for k=Nt-2:-1:1
                flds = (z[:,k+1],u[:,k:k+1],γ[k+1],ζ)
                z[:,k] = RRK_update_adj(flds,ops,dt,rk)
            end
            out_flds = z
        end
    else
    #RRK algorithm
        Nt_tmp = Nt_tmp + ceil(Int,Nt_tmp)
        f,η,∇η = ops
        u0 = in_flds

        u = zeros(length(u0),Nt_tmp)
        u[:,1] = u0
        γ = ones(Nt_tmp)

        t  = zeros(Nt_tmp)
        t[1] = t0
        k  = 0
        tk = t[1]
        dt_corr = dt #corrected time step size
        last_step = false

        while tk<T && abs(T-tk)/T>1e-13 && ~last_step
            k += 1
            if t[k]+dt > T
                # @show k
                # @show tk
                # @show t[k]+dt
                last_step = true
                dt_corr = T-t[k]
            end
            flds = (u[:,k],γ[k])

            u[:,k+1],γ[k+1] = RRK_update(flds,(f,η,∇η),dt_corr,rk)
            t[k+1] = t[k] + γ[k+1]*dt_corr
            tk = t[k+1]
        end

        Nt  = k+1
        u   = u[:,1:Nt]
        γ   = γ[1:Nt]

        out_flds = (u,γ)

        if return_Δη
            Δηu = zeros(Nt)
            ηu0 = η(u[:,1])
            for k=1:Nt-1
                Δηu[k+1] = η(u[:,k+1])-ηu0
            end
            out_flds = (out_flds...,Δηu)
        end
        if return_time
            t  = t[1:Nt]
            out_flds = (out_flds...,t,dt_corr)
        end
    end

    return out_flds
end



# Interpolation code
# function RK_interp_lin(x_in,y_in,x_out)
#     N = 1
#     if size(y_in,2)>1
#         N = size(y_in,1)
#     end
#     K = length(x_out)
#     y_out = zeros(N,K)
#     nodes = (x_in,)
#
#     if N==1
#         itp = interpolate(nodes,y_in,Gridded(Linear()))
#         for k=1:K-1
#             y_out[k] = itp(x_out[k])
#         end
#         y_out[K] = y_in[K]
#     else
#         for i=1:N
#             itp = interpolate(nodes,y_in[i,:],Gridded(Linear()))
#             for k=1:K-1
#                 y_out[i,k] = itp(x_out[k])
#             end
#             y_out[i,K] = y_in[i,K]
#         end
#     end
#     return y_out
# end
