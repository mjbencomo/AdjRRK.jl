## Algorithms for diagonally implicit RRK, assuming f(y) in y'=f(y) is linear,
#  i.e, f(y) = f*y where f is a matrix.

function DIRK_update(u,f,dt,rk::RK_struct)
    #assuming f is a matrix
    @unpack b,A = rk
    s = rk.stages
    U = zeros(length(u),s)
    Imat = Matrix{Float64}(I,length(u),length(u))

    U[:,1] = u
    U[:,1] = (Imat - dt*A[1,1]*f)\U[:,1]
    du = b[1].*f*U[:,1]
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f*U[:,j]
        end
        U[:,i] = u + dt*U[:,i]
        U[:,i] = ( Imat - dt*A[i,i]*f )\U[:,i]
        du += b[i]*f*U[:,i]
    end
    return @. u + dt*du
end

function DIRK_update_lin(flds,ops,dt,rk::RK_struct)
    #assuming f is a matrix
    @unpack b,A = rk
    s = rk.stages
    w,u  = flds
    f,df = ops
    Imat = Matrix{Float64}(I,length(u),length(u))

    #recomputing internal fwd stages
    U = zeros(length(u),s)
    U[:,1] = u
    U[:,1] = ( Imat - dt*A[1,1]*f )\U[:,1]
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f*U[:,j]
        end
        U[:,i] = u + dt*U[:,i]
        U[:,i] = ( Imat - dt*A[i,i]*f )\U[:,i]
    end

    #linear internal stages
    W = zeros(length(u),s)
    W[:,1] = w
    W[:,1] = ( Imat - dt*A[1,1]*f )\W[:,1]
    dw = b[1].*df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*df(U[:,j],W[:,j])
        end
        W[:,i] = w + dt*W[:,i]
        W[:,i] = ( Imat - dt*A[i,i]*f )\W[:,i]
        dw += b[i]*df(U[:,i],W[:,i])
    end
    return @. w + dt*dw
end

function DIRK_update_adj(flds,ops,dt,rk::RK_struct)
    #assuming problem is linear, f=matrix
    @unpack b,A = rk
    s = rk.stages
    z,u  = flds
    f,df = ops
    fT = transpose(f)
    Imat = Matrix{Float64}(I,length(u),length(u))

    #recomputing fwd internal stages
    U = zeros(length(u),s)
    U[:,1] = u
    U[:,1] = (Imat - dt*A[1,1]*f)\U[:,1]
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*f*U[:,j]
        end
        U[:,i] = u + dt*U[:,i]
        U[:,i] = ( Imat - dt*A[i,i]*f )\U[:,i]
    end

    #adjoint internal stages
    Z = zeros(length(u),s)
    Z[:,s] = b[s]*z
    Z[:,s] = dt*df(U[:,s],Z[:,s];adj=true)
    Z[:,s] = ( Imat - dt*A[s,s]*fT )\Z[:,s]
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = b[i]*z + Z[:,i]
        Z[:,i] = dt*df(U[:,i],Z[:,i];adj=true)
        Z[:,i] = ( Imat - dt*A[i,i]*fT )\Z[:,i]
        dz += Z[:,i]
    end
    return @. z + dz
end

function DIRK_fwd!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack t0,T,dt = ts
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)
    @pack! ts = dt,Nt

    @unpack f = arrks
    @unpack u0 = arrks

    u = zeros(length(u0),Nt)
    u[:,1] = u0

    for k=1:Nt-1
        u[:,k+1] = DIRK_update(u[:,k],f,dt,rk)
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

function DIRK_lin!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack u0_lin,u = arrks

    u_lin = zeros(size(u))
    u_lin[:,1] = u0_lin

    for k=1:Nt-1
        flds = (u_lin[:,k],u[:,k])
        u_lin[:,k+1] = DIRK_update_lin(flds,(f,df),dt,rk)
    end
    @pack! arrks = u_lin
end

function DIRK_adj!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack uT_adj,u = arrks

    u_adj = zeros(size(u))
    u_adj[:,Nt] = uT_adj

    for k=Nt-1:-1:1
        flds = (u_adj[:,k+1],u[:,k])
        u_adj[:,k] = DIRK_update_adj(flds,(f,df),dt,rk)
    end
    @pack! arrks = u_adj
end

function DIRK_solver!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct;lin=false,adj=false)
    if lin && ~(adj)
        DIRK_lin!(arrks,ts,rk)
    elseif adj
        DIRK_adj!(arrks,ts,rk)
    else
        DIRK_fwd!(arrks,ts,rk)
    end
end

function DIRK_int(u,f,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    U = zeros(length(u),s)
    F = zeros(length(u),s)
    Imat = Matrix{Float64}(I,length(u),length(u))

    U[:,1] = (Imat - dt*A[1,1]*f)\u
    F[:,1] = f*U[:,1]
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*F[:,j]
        end
        U[:,i] =  u + dt*U[:,i]
        U[:,i] = ( Imat - dt*A[i,i]*f )\U[:,i]
        F[:,i] = f*U[:,i]
    end
    return U,F
end

function DIRK_int_lin(flds,ops,dt,rk::RK_struct)
    #assuming linear problem
    @unpack b,A = rk
    s = rk.stages
    w,u  = flds
    f,df = ops
    Imat = Matrix{Float64}(I,length(u),length(u))
    U = zeros(length(u),s)
    F = zeros(length(u),s)
    W = zeros(length(u),s)
    JW = zeros(length(u),s)

    #computing internal fwd stages
    U,F = DIRK_int(u,f,dt,rk)

    #linear internal stages
    W[:,1] = w
    W[:,1] = ( Imat - dt*A[1,1]*f )\W[:,1]
    JW[:,1] = df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*JW[:,j]
        end
        W[:,i] = w + dt*W[:,i]
        W[:,i] = ( Imat - dt*A[i,i]*f )\W[:,i]
        JW[:,i] = df(U[:,i],W[:,i])
    end
    return U,F,W,JW
end

## DIRRK functions
function DIRRKγ0_update_lin(flds,ops,dt,rk::RK_struct)
    w,u,γ = flds
    new_b = γ .* rk.b
    new_rk = RK_struct(rk.stages,new_b,rk.A)
    return DIRK_update_lin((w,u),ops,dt,new_rk)
end

function DIRRKγ0_update_adj(flds,ops,dt,rk::RK_struct)
    z,u,γ = flds
    new_b = γ .* rk.b
    new_rk = RK_struct(rk.stages,new_b,rk.A)
    return DIRK_update_adj((z,u),ops,dt,new_rk)
end

function DIRRK_update(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    u,γ0   = flds
    f,η,∇η = ops
    U,F = DIRK_int(u,f,dt,rk)

    d = zeros(length(u))
    e = 0
    for j=1:s
        d += b[j].*F[:,j]
        e += b[j]*( ∇η(U[:,j])⋅F[:,j] )
    end
    d *= dt
    e *= dt
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
        imax = 10
        while r(γL)>0 && i<=imax
            γR = γL
            γL = γL*0.5
            i  = i+1
            if(i==imax)
                @show dt
                @show γ0
                @show γL
                @show r(γ0)
                @show r(γL)
                @show r(0)
                @show r(0.5)
                @show r(1)
                error("Exiting from RRK_update, could not find γL<0 under 10 iterations")
            end
        end

        i = 1
        while r(γR)<0 && i<=imax
            γL = γR
            γR = γR*2
            i  = i+1
            if i==imax
                error("Exiting from RRK_update, could not find γR>0 under 10 iterations")
            end
        end
        γ = AdjRRK.bisect(γL,γR,r)
    end
    u = u + γ*d

    return u,γ
end

function DIRRK_update_lin(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u,γ,ϱ = flds
    f,df,∇η,Hη = ops

    u1 = u[:,1] #this may need to be changed
    u2 = u[:,2]

    U,F,W,JW = DIRK_int_lin((w,u1),(f,df),dt,rk)

    ∇ηu2 = ∇η(u2)
    ∇γ_uw = ( ∇ηu2 - ∇η(u1) )⋅w
    ∇γ_UW = 0
    ∂r_γ  = 0
    for j=1:s
        ∇ηd = ∇ηu2-∇η(U[:,j])
        ∂r_γ += b[j]*( ∇ηd⋅F[:,j] )
        ∇γ_UW += b[j]*( ∇ηd⋅JW[:,j] - F[:,j]⋅Hη(U[:,j],W[:,j]) )
    end
    ∂r_γ *= dt
    ∇γ_UW *= γ*dt
    ρ = -(∇γ_UW + ∇γ_uw)/∂r_γ

    dw = zeros(length(w))
    for j=1:s
        dw += b[j]*( γ*JW[:,j] + ρ*F[:,j] )
    end
    ϱ += ρ
    return (@. w + dt*dw, ϱ )
end

function DIRRK_update_adj(flds,ops,dt,rk::RK_struct)
    #assumes f=matrix
    @unpack b,A = rk
    s = rk.stages
    z,u,γ,ζ = flds
    f,df,∇η,Hη = ops
    fT = transpose(f)
    Imat = Matrix{Float64}(I,length(z),length(z))

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = DIRK_int(u1,f,dt,rk)

    ∇ηu2 = ∇η(u2)
    ∇γ_U = zeros(length(u1),s)
    ∂r_γ = 0
    ξ  = 0
    for j=1:s
        ∇ηd       = ∇ηu2-∇η(U[:,j])
        ∇γ_U[:,j] = b[j].*( df(U[:,j],∇ηd;adj=true) - Hη(U[:,j],F[:,j];adj=true) )
        ∂r_γ     += b[j]*( ∇ηd⋅F[:,j] )
        ξ        += b[j]*( F[:,j]⋅z )
    end
    ∂r_γ *= dt
    ∇γ_U *= -γ*dt/∂r_γ
    ∇γ_u  = ( ∇η(u1)-∇ηu2 )./∂r_γ
    ξ    *= dt

    #adjoint internal stages
    Z = zeros(length(z),s)
    Z[:,s] = γ*b[s]*z
    Z[:,s] = dt*df(U[:,s],Z[:,s];adj=true) + (ξ-ζ)*∇γ_U[:,s]
    Z[:,s] = ( Imat - dt*A[s,s]*fT )\Z[:,s]
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = γ*b[i]*z + Z[:,i]
        Z[:,i] = dt*df(U[:,i],Z[:,i];adj=true) + (ξ-ζ)*∇γ_U[:,i]
        Z[:,i] = ( Imat - dt*A[i,i]*fT )\Z[:,i]
        dz += Z[:,i]
    end
    return @. z + dz + (ξ-ζ)*∇γ_u
end

function DIRRK_update_lin_last(flds,ops,dt_corr,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u,γ,ϱ,dt_old = flds
    f,df,∇η,Hη = ops
    Imat = Matrix{Float64}(I,length(w),length(w))

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = DIRK_int(u1,f,dt_corr,rk)

    W = zeros(length(w),s)
    JW = zeros(length(w),s)
    W[:,1] = ( Imat - dt_corr*A[1,1]*f )\w
    JW[:,1] = df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*( JW[:,j] - (dt_old/dt_corr*ϱ)*F[:,j] )
        end
        W[:,i] = w + dt_corr*W[:,i] - (dt_old*ϱ*A[i,i])*F[:,i]
        W[:,i] = ( Imat - (dt_corr*A[i,i])*f )\W[:,i]
        JW[:,i] = df(U[:,i],W[:,i])
    end

    ∇ηu2 = ∇η(u2)
    ∇γ_uw = ( ∇ηu2 - ∇η(u1) )⋅w
    ∇γ_UW = 0
    ∂r_γ  = 0
    for i=1:s
        ∇ηd    = ∇ηu2-∇η(U[:,i])
        ∂r_γ  += b[i]*( ∇ηd⋅F[:,i] )
        ∇γ_UW += b[i]*( ∇ηd⋅JW[:,i] - F[:,i]⋅Hη(U[:,i],W[:,i]) )
    end
    ∂r_γ *= dt_corr
    ∇γ_UW *= γ*dt_corr
    ρ = -( ∇γ_UW + ∇γ_uw )/∂r_γ

    dw = zeros(length(w))
    for i=1:s
        dw += b[i]*( γ*JW[:,i] + ρ*F[:,i] )
    end
    return @. w + dt_corr*dw
end

function DIRRK_update_adj_last(flds,ops,dt_corr,rk::RK_struct)
    #assuming f=matrix
    @unpack b,A = rk
    s = rk.stages
    z,u,γ,dt_old = flds
    f,df,∇η,Hη = ops
    fT = transpose(f)
    Imat = Matrix{Float64}(I,length(z),length(z))

    u1 = u[:,1]
    u2 = u[:,2]

    #recomputing internal fwd stages
    U,F = DIRK_int(u1,f,dt_corr,rk)

    ∇ηu2 = ∇η(u2)
    ∇γ_U = zeros(length(u1),s)
    ∂r_γ = 0
    ξ    = 0
    for j=1:s
        ∇ηd       = ∇ηu2-∇η(U[:,j])
        ∇γ_U[:,j] = b[j].*( df(U[:,j],∇ηd;adj=true) - Hη(U[:,j],F[:,j];adj=true) )
        ∂r_γ     += b[j]*( ∇ηd⋅F[:,j] )
        ξ        += b[j]*( F[:,j]⋅z )
    end
    ∂r_γ *= dt_corr
    ∇γ_U *= -γ*dt_corr/∂r_γ
    ∇γ_u  = ( ∇η(u1)-∇ηu2 )./∂r_γ
    ξ    *= dt_corr

    #adjoint internal stages
    Z = zeros(length(z),s)
    Z[:,s] = γ*b[s]*z
    Z[:,s] = dt_corr*df(U[:,s],Z[:,s];adj=true) + ξ*∇γ_U[:,s]
    Z[:,s] = ( Imat - dt_corr*A[s,s]*fT )\Z[:,s]
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = γ*b[i]*z + Z[:,i]
        Z[:,i] = dt_corr*df(U[:,i],Z[:,i];adj=true) + ξ*∇γ_U[:,i]
        Z[:,i] = ( Imat - dt_corr*A[i,i]*fT )\Z[:,i]
        dz += Z[:,i]
    end

    ζ  = 0
    for j=1:s
        for i=1:j
            ζ += A[j,i]*( F[:,i]⋅Z[:,j] )
        end
    end
    ζ *= dt_old

    return (@. z + dz + ξ*∇γ_u, ζ)
end

function DIRRK_fwd!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack t0,T,dt = ts
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)
    @pack! ts = dt
    Nt_tmp = Nt + ceil(Int,Nt/2)

    @unpack f,η,∇η = arrks
    @unpack u0 = arrks

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

        if(k+1>Nt_tmp)
            err("Temp data size not large enough. dt must be not small enough, resulting in |γ-1|>>0")
        end

        u[:,k+1],γ[k+1] = DIRRK_update(flds,(f,η,∇η),dt_corr,rk)
        t[k+1] = t[k] + γ[k+1]*dt_corr
        tk = t[k+1]
    end

    Nt  = k+1
    u   = u[:,1:Nt]
    γ   = γ[1:Nt]

    @pack! arrks = u,γ
    @pack! ts = Nt,dt_corr

    if arrks.return_time
        t  = t[1:Nt]
        @pack! ts = t
    end
    if arrks.return_Δη
        Δη = zeros(Nt)
        η0 = η(u[:,1])
        for k=1:Nt-1
            Δη[k+1] = η(u[:,k+1])-η0
        end
        @pack! arrks = Δη
    end
end

function DIRRK_lin!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt,dt_corr = ts
    @unpack f,df = arrks
    @unpack u0_lin,u,γ = arrks

    u_lin = zeros(size(u))
    u_lin[:,1] = u0_lin

    if arrks.γ_cnst
        for k=1:Nt-2
            flds = (u_lin[:,k],u[:,k],γ[k+1])
            u_lin[:,k+1] = DIRRKγ0_update_lin(flds,(f,df),dt,rk)
        end
        #last step
        flds = (u_lin[:,Nt-1],u[:,Nt-1],γ[Nt])
        u_lin[:,Nt] = DIRRKγ0_update_lin(flds,(f,df),dt_corr,rk)
    else
        @unpack ∇η,Hη = arrks
        ϱ = 0
        for k=1:Nt-2
            flds = (u_lin[:,k],u[:,k:k+1],γ[k+1],ϱ)
            u_lin[:,k+1],ϱ = DIRRK_update_lin(flds,(f,df,∇η,Hη),dt,rk)
        end

        if arrks.dt_cnst
             ϱ = 0
        end
        flds = (u_lin[:,Nt-1],u[:,Nt-1:Nt],γ[Nt],ϱ,dt)
        u_lin[:,Nt] = DIRRK_update_lin_last(flds,(f,df,∇η,Hη),dt_corr,rk)
    end
    @pack! arrks = u_lin
end

function DIRRK_adj!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct)
    @unpack dt,Nt,dt_corr = ts
    @unpack f,df = arrks
    @unpack uT_adj,u,γ = arrks

    u_adj = zeros(size(u))
    u_adj[:,Nt] = uT_adj

    if arrks.γ_cnst
        flds = (u_adj[:,Nt],u[:,Nt-1],γ[Nt])
        u_adj[:,Nt-1] = DIRRKγ0_update_adj(flds,(f,df),dt_corr,rk)
        for k=Nt-2:-1:1
            flds = (u_adj[:,k+1],u[:,k],γ[k+1])
            u_adj[:,k] = DIRRKγ0_update_adj(flds,(f,df),dt,rk)
        end
    else
        @unpack ∇η,Hη = arrks
        flds = (u_adj[:,Nt],u[:,Nt-1:Nt],γ[Nt],dt)
        u_adj[:,Nt-1],ζ = DIRRK_update_adj_last(flds,(f,df,∇η,Hη),dt_corr,rk)

        if arrks.dt_cnst
            ζ = 0
        end
        for k=Nt-2:-1:1
            flds = (u_adj[:,k+1],u[:,k:k+1],γ[k+1],ζ)
            u_adj[:,k] = DIRRK_update_adj(flds,(f,df,∇η,Hη),dt,rk)
        end
    end
    @pack! arrks = u_adj
end

function DIRRK_solver!(arrks::AdjRRK_struct,ts::Time_struct,rk::RK_struct;lin=false,adj=false)
    if lin && ~(adj)
        DIRRK_lin!(arrks,ts,rk)
    elseif adj
        DIRRK_adj!(arrks,ts,rk)
    else
        DIRRK_fwd!(arrks,ts,rk)
    end
end
