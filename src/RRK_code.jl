## Algorithms for simple explicit RRK methods.
# function RK_int(u,f,dt,rk::RKs)
#     S = rk.stages
#     U = zeros(length(u),S)
#     F = zeros(length(u),S)
#
#     U[:,1] = u
#     F[:,1] = f(u)
#     for s=2:S
#         U[:,s] = @. u + dt*rk.a[s-1]*F[:,s-1]
#         F[:,s] = f(U[:,s])
#     end
#     return U,F
# end
#
# function RK_int_lin(flds,ops,dt,rk::RKs)
#     w,u  = flds
#     f,df = ops
#     S = rk.stages
#     U = zeros(length(u),S)
#     F = zeros(length(u),S)
#     W = zeros(length(u),S)
#     JW = zeros(length(u),S)
#
#     U[:,1] = u
#     F[:,1] = f(U[:,1])
#     W[:,1] = w
#     JW[:,1] = df(U[:,1],W[:,1])
#     for s=2:S
#         U[:,s] = @. u + dt*rk.a[s-1]*F[:,s-1]
#         F[:,s] = f(U[:,s])
#
#         W[:,s] = @. w + dt*rk.a[s-1]*JW[:,s-1]
#         JW[:,s] = df(U[:,s],W[:,s])
#     end
#     return U,F,W,JW
# end
#
# function RRK_update(flds,ops,dt,rk::RKs)
#     u,γ0   = flds
#     f,η,∇η = ops
#     U,F = RK_int(u,f,dt,rk)
#
#     d = zeros(length(u))
#     e = 0
#     for s=1:rk.stages
#         d = @. d + rk.b[s].*F[:,s]
#         e = e + rk.b[s]*( ∇η(U[:,s])⋅F[:,s] )
#     end
#     d = @. dt*d
#     e = dt*e
#     ηu = η(u)
#
#     function r(γ)
#         return η(@. u + γ*d) - ηu - γ*e
#     end
#
#     γ = 1
#     #this tol should be relative to η(u0)
#     if abs(r(1)) > 1e-15
#         γL = 0.9*γ0
#         γR = 1.1*γ0
#         i = 1
#         imax = 10
#         while r(γL)>0 && i<=imax
#             γR = γL
#             γL = γL*0.5
#             i  = i+1
#             if(i==imax)
#                 @show dt
#                 @show γ0
#                 @show γL
#                 @show r(γ0)
#                 @show r(γL)
#                 @show r(0)
#                 @show r(0.5)
#                 @show r(1)
#                 error("Exiting from RRK_update, could not find γL<0 under 10 iterations")
#             end
#         end
#
#         i = 1
#         while r(γR)<0 && i<=imax
#             γL = γR
#             γR = γR*2
#             i  = i+1
#             if i==imax
#                 error("Exiting from RRK_update, could not find γR>0 under 10 iterations")
#             end
#         end
#         γ = bisect(γL,γR,r)
#     end
#     u = u + γ*d
#
#     return u,γ
# end
#
# function RRK_update_lin(flds,ops,dt,rk::RKs)
#     w,u,γ,ϱ = flds
#     f,df,∇η,Hη = ops
#     S = rk.stages
#     b = rk.b
#     a = rk.a
#
#     u1 = u[:,1]
#     u2 = u[:,2]
#
#     U,F,W,JW = RK_int_lin((w,u1),(f,df),dt,rk)
#
#     ∇ηu2 = ∇η(u2)
#     Γw = ( ∇ηu2 - ∇η(u1) )⋅w
#     Γ̃W = 0
#     rγ = 0
#     for s=1:S
#         ∇ηd = ∇ηu2-∇η(U[:,s])
#         rγ = rγ + b[s]*( ∇ηd⋅F[:,s] )
#         Γ̃W = Γ̃W + b[s]*( ∇ηd⋅JW[:,s] - F[:,s]⋅Hη(U[:,s],W[:,s]) )
#     end
#     rγ = dt*rγ
#     Γ̃W = γ*dt*Γ̃W
#     ρ = Γ̃W + Γw
#     ρ = -ρ/rγ
#
#     dw = zeros(length(w))
#     for s=1:S
#         dw = @. dw + b[s]*( γ*JW[:,s] + ρ*F[:,s] )
#     end
#     ϱ += ρ
#     return (@. w + dt*dw, ϱ )
# end
#
# function RRK_update_adj(flds,ops,dt,rk::RKs)
#     z,u,γ,ζ = flds
#     f,df,∇η,Hη = ops
#     S = rk.stages
#     b = rk.b
#     a = rk.a
#
#     u1 = u[:,1]
#     u2 = u[:,2]
#
#     U,F = RK_int(u1,f,dt,rk)
#
#     ∇ηu2 = ∇η(u2)
#     Γ̃  = zeros(length(u1),S)
#     rγ = 0
#     ξ  = 0
#     for s=1:S
#         ∇ηd    = ∇ηu2-∇η(U[:,s])
#         Γ̃[:,s] = b[s].*( df(U[:,s],∇ηd;adj=true) - Hη(U[:,s],F[:,s];adj=true) )
#         rγ     = rγ + b[s]*( ∇ηd⋅F[:,s] )
#         ξ      = ξ  + b[s]*( F[:,s]⋅z )
#     end
#     rγ = dt*rγ
#     Γ̃  = -γ*dt/rγ .* Γ̃
#     Γ  = ( ∇η(u1)-∇ηu2 )./rγ
#     ξ  = dt*ξ - ζ
#
#     Z = df(U[:,S],z;adj=true)
#     Z = @. dt*γ*b[S]*Z + ξ*Γ̃[:,S]
#     dz = Z
#     for s=S-1:-1:1
#         Z = @. γ*b[s]*z + a[s]*Z
#         Z = df(U[:,s],Z;adj=true)
#         Z = @. dt*Z + ξ*Γ̃[:,s]
#         dz = @. dz + Z
#     end
#     return @. z + dz + ξ*Γ
# end
#
# function RRK_update_lin_last(flds,ops,dt_corr,rk::RKs)
#     w,u,γ,ϱ,dt_old = flds
#     f,df,∇η,Hη = ops
#     S = rk.stages
#     b = rk.b
#     a = rk.a
#
#     u1 = u[:,1]
#     u2 = u[:,2]
#
#     U = zeros(length(u1),S)
#     F = zeros(length(u1),S)
#     W = zeros(length(u1),S)
#     JW = zeros(length(u1),S)
#
#     U[:,1] = u1
#     F[:,1] = f(U[:,1])
#     W[:,1] = w
#     JW[:,1] = df(U[:,1],W[:,1])
#     for s=2:S
#         U[:,s] = @. u1 + dt_corr*a[s-1]*F[:,s-1]
#         F[:,s] = f(U[:,s])
#         W[:,s] = @. w + dt_corr*a[s-1]*( JW[:,s-1] - dt_old/dt_corr*ϱ*F[:,s-1] )
#         JW[:,s] = df(U[:,s],W[:,s])
#     end
#
#     ∇ηu2 = ∇η(u2)
#     Γw = ( ∇ηu2 - ∇η(u1) )⋅w
#     Γ̃W = 0
#     rγ = 0
#     for s=1:S
#         ∇ηd = ∇ηu2-∇η(U[:,s])
#         rγ = rγ + b[s]*( ∇ηd⋅F[:,s] )
#         Γ̃W = Γ̃W + b[s]*( ∇ηd⋅JW[:,s] - F[:,s]⋅Hη(U[:,s],W[:,s]) )
#     end
#     rγ = dt_corr*rγ
#     Γ̃W = γ*dt_corr*Γ̃W
#     ρ = Γ̃W + Γw
#     ρ = -ρ/rγ
#
#     dw = zeros(length(w))
#     for s=1:S
#         dw = @. dw + b[s]*( γ*JW[:,s] + ρ*F[:,s] )
#     end
#     return @. w + dt_corr*dw
# end
#
# function RRK_update_adj_last(flds,ops,dt_corr,rk::RKs)
#     z,u,γ,dt_old = flds
#     f,df,∇η,Hη = ops
#     S = rk.stages
#     b = rk.b
#     a = rk.a
#
#     u1 = u[:,1]
#     u2 = u[:,2]
#
#     U,F = RK_int(u1,f,dt_corr,rk)
#
#     ∇ηu2 = ∇η(u2)
#     Γ̃  = zeros(length(u1),S)
#     rγ = 0
#     ξ  = 0
#     for s=1:S
#         ∇ηd    = ∇ηu2-∇η(U[:,s])
#         Γ̃[:,s] = b[s].*( df(U[:,s],∇ηd;adj=true) - Hη(U[:,s],F[:,s];adj=true) )
#         rγ     = rγ + b[s]*( ∇ηd⋅F[:,s] )
#         ξ      = ξ  + b[s]*( F[:,s]⋅z )
#     end
#     rγ = dt_corr*rγ
#     Γ̃  = -γ*dt_corr/rγ .* Γ̃
#     Γ  = ( ∇η(u1)-∇ηu2 )./rγ
#     ξ  = dt_corr*ξ
#
#     Z = df(U[:,S],z;adj=true)
#     Z = @. dt_corr*γ*b[S]*Z + ξ*Γ̃[:,S]
#     dz = Z
#     ζ = 0
#     for s=S-1:-1:1
#         ζ += a[s]*( F[:,s]⋅Z )
#         Z = @. γ*b[s]*z + a[s]*Z
#         Z = df(U[:,s],Z;adj=true)
#         Z = @. dt_corr*Z + ξ*Γ̃[:,s]
#         dz = @. dz + Z
#     end
#     ζ *= dt_old
#     return (@. z + dz + ξ*Γ, ζ)
# end



function RK_int(u,f,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    U = zeros(length(u),s)
    F = zeros(length(u),s)

    U[:,1] = u
    F[:,1] = f(u)
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*F[:,j]
        end
        U[:,i] = u + dt*U[:,i]
        F[:,i] = f(U[:,i])
    end
    return U,F
end

function RK_int_lin(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u  = flds
    f,df = ops
    U = zeros(length(u),s)
    F = zeros(length(u),s)
    W = zeros(length(u),s)
    JW = zeros(length(u),s)

    #computing internal fwd stages
    U[:,1] = u
    F[:,1] = f(U[:,1])
    for i=2:s
        for j=1:i-1
            U[:,i] += A[i,j]*F[:,j]
        end
        U[:,i] = u + dt*U[:,i]
        F[:,i] = f(U[:,i])
    end

    #linear internal stages
    W[:,1] = w
    JW[:,1] = df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*JW[:,j]
        end
        W[:,i] = w + dt*W[:,i]
        JW[:,i] = df(U[:,i],W[:,i])
    end
    return U,F,W,JW
end

function bisect(γL,γR,r;rtol=1e-12,γtol=1e-12,nmax=10000)
    γm = 0
    rm = 0
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

function RRK_update(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    u,γ0   = flds
    f,η,∇η = ops
    U,F = RK_int(u,f,dt,rk)

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
        γ = bisect(γL,γR,r)
    end
    u = u + γ*d

    return u,γ
end

function RRKγ0_update_lin(flds,ops,dt,rk::RK_struct)
    w,u,γ = flds
    new_b = γ .* rk.b
    new_rk = RK_struct(rk.stages,new_b,rk.A)
    return RK_update_lin((w,u),ops,dt,new_rk)
end

function RRKγ0_update_adj(flds,ops,dt,rk::RK_struct)
    z,u,γ = flds
    new_b = γ .* rk.b
    new_rk = RK_struct(rk.stages,new_b,rk.A)
    return RK_update_adj((z,u),ops,dt,new_rk)
end

function RRK_update_lin(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u,γ,ϱ = flds
    f,df,∇η,Hη = ops

    u1 = u[:,1] #this may need to be changed
    u2 = u[:,2]

    U,F,W,JW = RK_int_lin((w,u1),(f,df),dt,rk)

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

function RRK_update_adj(flds,ops,dt,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    z,u,γ,ζ = flds
    f,df,∇η,Hη = ops

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = RK_int(u1,f,dt,rk)

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
    Z[:,s] = dt*γ*b[s]*df(U[:,s],z;adj=true) + (ξ-ζ)*∇γ_U[:,s]
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = γ*b[i]*z + Z[:,i]
        Z[:,i] = df(U[:,i],Z[:,i];adj=true)
        Z[:,i] = dt*Z[:,i] + (ξ-ζ)*∇γ_U[:,i]
        dz += Z[:,i]
    end
    return @. z + dz + (ξ-ζ)*∇γ_u
end

function RRK_update_lin_last(flds,ops,dt_corr,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    w,u,γ,ϱ,dt_old = flds
    f,df,∇η,Hη = ops

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = RK_int(u1,f,dt_corr,rk)

    W = zeros(length(w),s)
    JW = zeros(length(w),s)
    W[:,1] = w
    JW[:,1] = df(U[:,1],W[:,1])
    for i=2:s
        for j=1:i-1
            W[:,i] += A[i,j]*( JW[:,j] - dt_old/dt_corr*ϱ*F[:,j] )
        end
        W[:,i] = @. w + dt_corr*W[:,i]
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

function RRK_update_adj_last(flds,ops,dt_corr,rk::RK_struct)
    @unpack b,A = rk
    s = rk.stages
    z,u,γ,dt_old = flds
    f,df,∇η,Hη = ops

    u1 = u[:,1]
    u2 = u[:,2]

    U,F = RK_int(u1,f,dt_corr,rk)

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

    Z = zeros(length(z),s)
    Z[:,s] = dt_corr*γ*b[s]*df(U[:,s],z;adj=true) + ξ*∇γ_U[:,s]
    dz = Z[:,s]
    for i=s-1:-1:1
        for j=i+1:s
            Z[:,i] += A[j,i]*Z[:,j]
        end
        Z[:,i] = γ*b[i]*z + Z[:,i]
        Z[:,i] = df(U[:,i],Z[:,i];adj=true)
        Z[:,i] = dt_corr*Z[:,i] + ξ*∇γ_U[:,i]
        dz += Z[:,i]
    end

    ζ  = 0
    for j=2:s
        for i=1:j-1
            ζ += A[j,i]*( F[:,i]⋅Z[:,j] )
        end
    end
    ζ *= dt_old

    return (@. z + dz + ξ*∇γ_u, ζ)
end

## IDT SOLVER
#
# INPUTS:
#   arrks::AdjRRK_struct
#       [IDT fwd case]
#           u0,     initial condition
#           f,      RHS function
#           η,      entropy function
#           ∇η,     gradient of entropy
#           return_time,
#           return_Δη,
#       [IDT lin case; if lin=true]
#           u,      solution of forward RK
#           γ,      relaxation parameter
#           u0_lin, initial condition for linearized solution
#           f,      RHS function
#           df,     Jacobian of RHS function
#           ∇η,     gradient of entropy (if γ_cnst=false)
#           Hη,     Hessian of entropy (if γ_cnst=false)
#           γ_cnst, flag, true if γ should be const. when linearizing
#       [IDT adj case; if adj=true]
#           u,      solution of forward RK
#           γ,      relaxation parameter
#           uT_adj, final condition for adjoint solution
#           f,      RHS function
#           df,     Jacobian of RHS function
#           ∇η,     gradient of entropy (if γ_cnst=false)
#           Hη,     Hessian of entropy (if γ_cnst=false)
#           γ_cnst, flag, true if γ should be const. when linearizing
#   ts::Time_struct
#       [IDT fwd case]
#           t0,     intial time
#           T,      final time
#           dt,     time step size
#       [IDT lin case; if lin=true]
#           dt,     time step size
#           Nt,     number of time steps
#       [IDT adj case; if adj=true]
#           dt,     time step size
#           Nt,     number of time steps
#   rk::RKs
#
# OUTPUTS:
#   arrks::AdjRRK_struct
#       [IDT fwd case]
#           u,      solution
#           γ,      relaxation parameter
#           Δη,     change in entropy (if return_Δη=true)
#       [IDT lin case; if lin=true]
#           u_lin,  linearized solution
#       [IDT adj case; if adj=true]
#           u_adj,  adjoint solution
#   ts::Time_struct
#       [IDT fwd case]
#           dt,     may get modified
#           Nt,     number of time steps
#           t,      time grid points

function IDT_solver!(arrks::AdjRRK_struct,ts::Time_struct,rk;lin=false,adj=false)
    if lin && ~(adj)
        IDT_lin!(arrks,ts,rk)
    elseif adj
        IDT_adj!(arrks,ts,rk)
    else
        IDT_fwd!(arrks,ts,rk)
    end
end

function IDT_fwd!(arrks::AdjRRK_struct,ts::Time_struct,rk)
    @unpack t0,T,dt = ts
    Nt = ceil(Int,(T-t0)/dt)+1
    dt = (T-t0)/(Nt-1)
    @pack! ts = dt,Nt

    @unpack f,η,∇η = arrks
    @unpack u0 = arrks

    u = zeros(length(u0),Nt)
    u[:,1] = u0
    γ = ones(Nt)

    for k=1:Nt-1
        flds = (u[:,k],γ[k])
        u[:,k+1],γ[k+1] = RRK_update(flds,(f,η,∇η),dt,rk)
    end
    @pack! arrks = u,γ

    if arrks.return_time
        t = range(t0,stop=T,length=Nt) |> collect
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

function IDT_lin!(arrks::AdjRRK_struct,ts::Time_struct,rk)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack u0_lin,u,γ = arrks

    u_lin = zeros(size(u))
    u_lin[:,1] = u0_lin

    if arrks.γ_cnst
        for k=1:Nt-1
            flds = (u_lin[:,k],u[:,k],γ[k+1])
            u_lin[:,k+1] = RRKγ0_update_lin(flds,(f,df),dt,rk)
        end
    else
        @unpack ∇η,Hη = arrks
        ϱ = 0
        for k=1:Nt-1
            flds = (u_lin[:,k],u[:,k:k+1],γ[k+1],ϱ)
            u_lin[:,k+1],ϱ = RRK_update_lin(flds,(f,df,∇η,Hη),dt,rk)
        end
    end
    @pack! arrks = u_lin
end

function IDT_adj!(arrks::AdjRRK_struct,ts::Time_struct,rk)
    @unpack dt,Nt = ts
    @unpack f,df = arrks
    @unpack uT_adj,u,γ = arrks

    u_adj = zeros(size(u))
    u_adj[:,Nt] = uT_adj

    if arrks.γ_cnst
        for k=Nt-1:-1:1
            flds = (u_adj[:,k+1],u[:,k],γ[k+1])
            u_adj[:,k] = RRKγ0_update_adj(flds,(f,df),dt,rk)
        end
    else
        @unpack ∇η,Hη = arrks
        for k=Nt-1:-1:1
            flds = (u_adj[:,k+1],u[:,k:k+1],γ[k+1],0)
            u_adj[:,k] = RRK_update_adj(flds,(f,df,∇η,Hη),dt,rk)
        end
    end
    @pack! arrks = u_adj
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
function RRK_solver!(arrks::AdjRRK_struct,ts::Time_struct,rk;lin=false,adj=false)
    if lin && ~(adj)
        RRK_lin!(arrks,ts,rk)
    elseif adj
        RRK_adj!(arrks,ts,rk)
    else
        RRK_fwd!(arrks,ts,rk)
    end
end

function RRK_fwd!(arrks::AdjRRK_struct,ts::Time_struct,rk)
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

        u[:,k+1],γ[k+1] = RRK_update(flds,(f,η,∇η),dt_corr,rk)
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

function RRK_lin!(arrks::AdjRRK_struct,ts::Time_struct,rk)
    @unpack dt,Nt,dt_corr = ts
    @unpack f,df = arrks
    @unpack u0_lin,u,γ = arrks

    u_lin = zeros(size(u))
    u_lin[:,1] = u0_lin

    if arrks.γ_cnst
        for k=1:Nt-2
            flds = (u_lin[:,k],u[:,k],γ[k+1])
            u_lin[:,k+1] = RRKγ0_update_lin(flds,(f,df),dt,rk)
        end
        #last step
        flds = (u_lin[:,Nt-1],u[:,Nt-1],γ[Nt])
        u_lin[:,Nt] = RRKγ0_update_lin(flds,(f,df),dt_corr,rk)
    else
        @unpack ∇η,Hη = arrks
        ϱ = 0
        for k=1:Nt-2
            flds = (u_lin[:,k],u[:,k:k+1],γ[k+1],ϱ)
            u_lin[:,k+1],ϱ = RRK_update_lin(flds,(f,df,∇η,Hη),dt,rk)
        end

        if arrks.dt_cnst
             ϱ = 0
        end
        flds = (u_lin[:,Nt-1],u[:,Nt-1:Nt],γ[Nt],ϱ,dt)
        u_lin[:,Nt] = RRK_update_lin_last(flds,(f,df,∇η,Hη),dt_corr,rk)
    end
    @pack! arrks = u_lin
end

function RRK_adj!(arrks::AdjRRK_struct,ts::Time_struct,rk)
    @unpack dt,Nt,dt_corr = ts
    @unpack f,df = arrks
    @unpack uT_adj,u,γ = arrks

    u_adj = zeros(size(u))
    u_adj[:,Nt] = uT_adj

    if arrks.γ_cnst
        flds = (u_adj[:,Nt],u[:,Nt-1],γ[Nt])
        u_adj[:,Nt-1] = RRKγ0_update_adj(flds,(f,df),dt_corr,rk)
        for k=Nt-2:-1:1
            flds = (u_adj[:,k+1],u[:,k],γ[k+1])
            u_adj[:,k] = RRKγ0_update_adj(flds,(f,df),dt,rk)
        end
    else
        @unpack ∇η,Hη = arrks
        flds = (u_adj[:,Nt],u[:,Nt-1:Nt],γ[Nt],dt)
        u_adj[:,Nt-1],ζ = RRK_update_adj_last(flds,(f,df,∇η,Hη),dt_corr,rk)

        if arrks.dt_cnst
            ζ = 0
        end
        for k=Nt-2:-1:1
            flds = (u_adj[:,k+1],u[:,k:k+1],γ[k+1],ζ)
            u_adj[:,k] = RRK_update_adj(flds,(f,df,∇η,Hη),dt,rk)
        end
    end
    @pack! arrks = u_adj
end

# function RRK_solver(in_flds,ops,Time,rk;
# lin=false,
# adj=false,
# γcnst=false,
# corr=true,
# return_time=false,
# return_Δη=false)
#
#     t0,T,dt = Time
#     Nt_tmp = ceil(Int,(T-t0)/dt)+1
#     dt = (T-t0)/(Nt_tmp-1)
#
#     if lin && ~(adj)
#     #linearized RRK algorithm
#         w0,u,γ,dt_corr = in_flds
#         Nt = length(γ)
#         w = zeros(size(u))
#         w[:,1] = w0
#
#         if γcnst
#             for k=1:Nt-2
#                 flds = (w[:,k],u[:,k],γ[k+1])
#                 w[:,k+1] = RRKγ0_update_lin(flds,ops,dt,rk)
#             end
#             #last step
#             flds = (w[:,Nt-1],u[:,Nt-1],γ[Nt])
#             w[:,Nt] = RRKγ0_update_lin(flds,ops,dt_corr,rk)
#             out_flds = w
#         else
#             ϱ = 0
#             for k=1:Nt-2
#                 flds = (w[:,k],u[:,k:k+1],γ[k+1],ϱ)
#                 w[:,k+1],ϱ = RRK_update_lin(flds,ops,dt,rk)
#             end
#
#             if ~corr
#                 ϱ = 0
#             end
#             flds = (w[:,Nt-1],u[:,Nt-1:Nt],γ[Nt],ϱ,dt)
#             w[:,Nt] = RRK_update_lin_last(flds,ops,dt_corr,rk)
#             out_flds = w
#         end
#
#     elseif adj
#     #adjoint RRK algorithm
#         zT,u,γ,dt_corr = in_flds
#         Nt = length(γ)
#         z = zeros(size(u))
#         z[:,Nt] = zT
#
#         if γcnst
#             flds = (z[:,Nt],u[:,Nt-1],γ[Nt])
#             z[:,Nt-1] = RRKγ0_update_adj(flds,ops,dt_corr,rk)
#             for k=Nt-2:-1:1
#                 flds = (z[:,k+1],u[:,k],γ[k+1])
#                 z[:,k] = RRKγ0_update_adj(flds,ops,dt,rk)
#             end
#             out_flds = z
#         else
#             flds = (z[:,Nt],u[:,Nt-1:Nt],γ[Nt],dt)
#             z[:,Nt-1],ζ = RRK_update_adj_last(flds,ops,dt_corr,rk)
#
#             if ~corr
#                 ζ = 0
#             end
#             for k=Nt-2:-1:1
#                 flds = (z[:,k+1],u[:,k:k+1],γ[k+1],ζ)
#                 z[:,k] = RRK_update_adj(flds,ops,dt,rk)
#             end
#             out_flds = z
#         end
#     else
#     #RRK algorithm
#         Nt_tmp = Nt_tmp + ceil(Int,Nt_tmp)
#         f,η,∇η = ops
#         u0 = in_flds
#
#         u = zeros(length(u0),Nt_tmp)
#         u[:,1] = u0
#         γ = ones(Nt_tmp)
#
#         t  = zeros(Nt_tmp)
#         t[1] = t0
#         k  = 0
#         tk = t[1]
#         dt_corr = dt #corrected time step size
#         last_step = false
#
#         while tk<T && abs(T-tk)/T>1e-13 && ~last_step
#             k += 1
#             if t[k]+dt > T
#                 # @show k
#                 # @show tk
#                 # @show t[k]+dt
#                 last_step = true
#                 dt_corr = T-t[k]
#             end
#             flds = (u[:,k],γ[k])
#
#             if(k+1>Nt_tmp)
#                 err("Temp data size not large enough. dt must be too big, resulting in |γ-1|>>0")
#             end
#
#             u[:,k+1],γ[k+1] = RRK_update(flds,(f,η,∇η),dt_corr,rk)
#             t[k+1] = t[k] + γ[k+1]*dt_corr
#             tk = t[k+1]
#         end
#
#         Nt  = k+1
#         u   = u[:,1:Nt]
#         γ   = γ[1:Nt]
#
#         out_flds = (u,γ)
#
#         if return_Δη
#             Δηu = zeros(Nt)
#             ηu0 = η(u[:,1])
#             for k=1:Nt-1
#                 Δηu[k+1] = η(u[:,k+1])-ηu0
#             end
#             out_flds = (out_flds...,Δηu)
#         end
#         if return_time
#             t  = t[1:Nt]
#             out_flds = (out_flds...,t,dt_corr)
#         end
#     end
#
#     return out_flds
# end



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
