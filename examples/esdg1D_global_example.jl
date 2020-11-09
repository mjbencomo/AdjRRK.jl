using Plots
using UnPack
using LinearAlgebra
using SparseArrays

using ForwardDiff
using StaticArrays

using StartUpDG
using StartUpDG.ExplicitTimestepUtils

using EntropyStableEuler.Fluxes1D
import EntropyStableEuler: γ
using FluxDiffUtils

include("HybridizedSBPUtils.jl")
using .HybridizedSBPUtils

using AdjRRK

N = 3
K1D = 16
CFL = .25
T = 0.5
Nfields = 3

################################################################################
## Initial DG stuff

rd = init_reference_interval(N)
VX,EToV = uniform_1D_mesh(K1D)
md = init_DG_mesh(VX,EToV,rd)

@unpack M,Dr,Vq,Pq,Vf,wf,nrJ = rd
Qr = Pq'*M*Dr*Pq
Ef = Vf*Pq
Br = diagm(wf.*nrJ)
Qrh = .5*[Qr-Qr' Ef'*Br;
        -Br*Ef  Br]
Vh = [Vq;Vf]
Ph = M\transpose(Vh)
VhP = Vh*Pq

# make skew symmetric versions of the operators"
Qrhskew = .5*(Qrh-transpose(Qrh))
SBP_ops = (Matrix(Qrhskew'),Vh,Ph,VhP)

LX = 2
build_periodic_boundary_maps!(md,rd,LX)

################################################################################
##

@unpack x = md
#initial condition
rho = @. 1 + .5*exp(-100*(x-.1)^2)
u = @. 0*x
p = @. rho^γ

Q0 = primitive_to_conservative(rho,u,p)
Q = primitive_to_conservative(rho,u,p)

# interpolate geofacs to both vol/surf nodes, pack back into md
@unpack rxJ = md
rxJ = [Vq;Vf]*rxJ # interp to hybridized points
@pack! md = rxJ

################################################################################
## build global matrices

Nh,Nq = size(VhP)
function rhs(uh,rd,md)
    @unpack rxJ,nxJ,mapP = md
    uP = ((x->x[mapP]) ∘ (x->x[Nq+1:Nh,:]))(uh)
    rhs = rxJ.*(Qrhskew*uh)
    rhs[Nq+1:Nh,:] += .5*uP.*nxJ
    return rhs
end

function rhsB(uh,rd,md)
    @unpack rxJ,nxJ,mapP = md
    uP = ((x->x[mapP]) ∘ (x->x[Nq+1:Nh,:]))(uh)
    rhs = zero.(rxJ.*(Qrhskew*uh))
    rhs[Nq+1:Nh,:] += .5*(uP.*nxJ)
    return rhs
end

Qg = droptol!(build_rhs_matrix(rhs,size(Vh,1),md.K,rd,md),1e-12)
Bg = droptol!(build_rhs_matrix(rhsB,size(Vh,1),md.K,rd,md),1e-12)
Vqg,Vhg,VhPg = (A->kron(I(md.K),A)).((Vq,Vh,VhP))
#Phg = kron(I(md.K),Ph)
Phg = kron(diagm(1 ./ md.J[1,:]),Ph)

QgTr = sparse(transpose(Qg))
global_SBP_ops = (Qg,QgTr,Bg,Vhg,Phg,VhPg,Vqg)

function LF_dissipation(rhoL,rhouL,EL,rhoR,rhouR,ER)
    QL = rhoL,rhouL,EL
    QR = rhoR,rhouR,ER
    λL = abs.(wavespeed_1D(QL...))
    λR = abs.(wavespeed_1D(QR...))
    λavg =.5*(λL+λR)
    # λavg = 1.
    return ((uL,uR)->.5*λavg*(uR-uL)).(QL,QR)
end

function rhs_global(Q,global_SBP_ops)
    Qx,QxTr,Bx,Vh,Ph,VhP,Vq = global_SBP_ops
    Qh = u_vfun((x->VhP*x).(v_ufun((x->Vq*x).(vec.(Q))...))...)
    Qhprim = conservative_to_primitive_beta(Qh...)
    rhsQ = hadamard_sum((QxTr,),tuple ∘ euler_fluxes, Qhprim) .-
            hadamard_sum((abs.(Bx),),tuple ∘ LF_dissipation, Qh)
    return (x->-2*Ph*x).(rhsQ)
end


## Jacobian stuff

# define Euler fluxes directly as functions of conservative variables
function F(UL,UR)
    # convert to flux variables
    function UtoQ(U)
        rho,rhou,E = U
        return (rho,rhou./rho,betafun(U...))
    end
    Fx = euler_fluxes(UtoQ(UL)...,UtoQ(UR)...)
    return SVector{length(Fx)}(Fx...)
end

# mappings between conservative and entropy variables and vice versa
dVdU_fun(U::SVector) = ForwardDiff.jacobian(U->SVector(v_ufun(U...)),U)
dUdV_fun(V::SVector) = ForwardDiff.jacobian(V->SVector(u_vfun(V...)),V)
# "slurping" individual arguments for compatibility with FluxDiffUtils.jl
dVdU_fun(U...) = dVdU_fun(SVector(U))
dUdV_fun(V...) = dUdV_fun(SVector(V))

# AD for jacobians
dFx(uL::SVector,uR::SVector) = ForwardDiff.jacobian(uR->F(uL,uR),uR)
dLF(uL::SVector,uR::SVector) = ForwardDiff.jacobian(
                uR->SVector(LF_dissipation(uL...,uR...)),uR)
# inputs are individual arguments, outputs wrapped in a tuple
dFx(rhoL,rhouL,EL,rhoR,rhouR,ER) = tuple(dFx(SVector(rhoL,rhouL,EL),SVector(rhoR,rhouR,ER)))
dLF(rhoL,rhouL,EL,rhoR,rhouR,ER) = tuple(dLF(SVector(rhoL,rhouL,EL),SVector(rhoR,rhouR,ER)))

# compute "parts" of Jacobians
function rhs_global_jac(Q,global_SBP_ops)

    Qx,QxTr,Bx,Vh,Ph,VhP,Vq = global_SBP_ops

    U = SVector(vec.(Q))
    Uq = (x->Vq*x).(U)
    VUh = SVector((x->VhP*x).(v_ufun(Uq...)))
    Qh = u_vfun(VUh...)
    dVdU = banded_matrix_function(dVdU_fun,Uq)
    dUdV = banded_matrix_function(dUdV_fun,VUh)
    Jblocks = hadamard_jacobian(tuple(Qx),:skew,dFx,Qh) # compute jacobian for entropy conservative part
    JLFblocks = hadamard_jacobian(tuple(abs.(Bx)),:skew,dLF,Qh) # compute jacobian for dissipation term

    # convert FluxDiffUtils.jl output (tuples of matrices) to global matrices
    dVdU,dUdV,Jac,JLF = (A->hcat(vcat.(A...)...)).((dVdU,dUdV,Jblocks,JLFblocks))

    # compute jacobian using chain rule
    return kron(I(Nfields),-2*Ph) * (Jac - JLF) * dUdV * kron(I(Nfields),VhP) * dVdU * kron(I(Nfields),Vq)
end



function tuple_to_arr(X)
    return vcat(X...)
end

function arr_to_tuple(X)
    Length = length(X)÷Nfields
    return ntuple(i->reshape(X,Length,Nfields)[:,i],Nfields)
end


################################################################################
##

function f(u)
    Q  = arr_to_tuple(u)
    fQ = rhs_global(Q,global_SBP_ops)
    return tuple_to_arr(fQ)
end

function df(u,δu;adj=false)
    Q = arr_to_tuple(u)
    J = rhs_global_jac(Q,global_SBP_ops)
    if adj
        return transpose(J)*δu
    end
    return J*δu
end

@unpack wJq = md
function η(u)
    Q = arr_to_tuple(u)
    Q = (x->reshape(x,length(rd.r),md.K)).(Q)
    return sum(wJq.*(Sfun((x->Vq*x).(Q)...)))
    #return sum(wJq.*(Vq*Sfun(Q...)))
end

function ∇η(u)
    Q = arr_to_tuple(u) #reshape into tuple of 1D arrays
    Q = (x->reshape(x,length(rd.r),md.K)).(Q) #reshape into tuple of 2D arrays

    #entropy projection step
    VU = v_ufun((x->Vq*x).(Q)...)
    VU = (x->Vq*Pq*x).(VU)
    VU = (x->wJq.*x).(VU)
    VU = (x->transpose(Vq)*x).(VU)
    VU = vec.(VU)
    return tuple_to_arr(VU)
end

function Hη(u,δu;adj=false)
    Q = arr_to_tuple(u)
    Qx,QxTr,Bx,Vh,Ph,VhP,Vq = global_SBP_ops

    WJq  = spdiagm(0 => wJq[:,1])
    WJqg = kron(I(md.K),WJq)
    Pqg  = kron(I(md.K),Pq)

    U = SVector(vec.(Q))
    Uq = (x->Vq*x).(U)
    dVdU = banded_matrix_function(dVdU_fun,Uq)

    dVdU = hcat(vcat.(dVdU...)...)

    H = kron(I(Nfields),transpose(Vq)*WJqg*Vq*Pqg) * dVdU * kron(I(Nfields),Vq)

    if adj
        return transpose(H)*δu
    end
    return H*δu
end

# #running derv test of df
# u = u0 #randn(size(u0)).+10
# δu = randn(size(u0))
# h0 = 2^(-10)
# Nref = 5
# h = zeros(Nref+1)
# h[1] = h0
# for n=1:Nref
#     h[n+1] = h[n]/2
# end
# errs = zeros(Nref+1)
# rate = zeros(Nref)
#
# fu = f(u)
# dfu = df(u,δu)
# for n=1:Nref+1
#     uh = u + h[n].*δu
#     fuh = f(uh)
#     dfuh = (fuh - fu) / h[n] #./norm(δu)
#     errs[n] = norm( dfuh - dfu )
# end
# @show errs
#
# for n=1:Nref
#     rate[n] = log2(errs[n]) - log2(errs[n+1])
# end
# @show rate

# #running derv test of ∇η
# u = u0#randn(size(u0)).+10
# δu = randn(size(u0))
# h0 = 2^(-8)
# Nref = 3
# h = zeros(Nref+1)
# h[1] = h0
# for n=1:Nref
#     h[n+1] = h[n]/2
# end
# errs = zeros(Nref+1)
# rate = zeros(Nref)
#
# ηu = η(u)
# ∇ηu = ∇η(u)⋅δu
# for n=1:Nref+1
#     uh = u + h[n].*δu
#     ηuh = η(uh)
#     ∇ηuh = (ηuh - ηu) / h[n] #./norm(δu)
#     errs[n] = norm( ∇ηuh - ∇ηu )
# end
# @show errs
#
# for n=1:Nref
#     rate[n] = log2(errs[n]) - log2(errs[n+1])
# end
# @show rate

# #running derv test of Hη
# u = randn(size(u0)).+10
# δu = randn(size(u0))
# h0 = 2^(-8)
# Nref = 3
# h = zeros(Nref+1)
# h[1] = h0
# for n=1:Nref
#     h[n+1] = h[n]/2
# end
# errs = zeros(Nref+1)
# rate = zeros(Nref)
#
# ∇ηu = ∇η(u)
# Hηu = Hη(u,δu)
# for n=1:Nref+1
#     uh = u + h[n].*δu
#     ∇ηuh = ∇η(uh)
#     Hηuh = (∇ηuh - ∇ηu) / h[n] #./norm(δu)
#     errs[n] = norm( Hηuh - Hηu )
# end
# @show errs
#
# for n=1:Nref
#     rate[n] = log2(errs[n]) - log2(errs[n+1])
# end
# @show rate

#time stuff
CN = (N+1)^2/2  # estimated trace constant
dt = CFL * 2 / (CN*K1D)
t0 = 0

Q0 = vec.(Q0)
u0 = tuple_to_arr(Q0)

# Q00 = (x->reshape(x,length(rd.r),md.K)).(Q0)
# display(scatter((x->rd.Vp*x).((x,Q00[3])),zcolor=rd.Vp*Q00[3],msw=0,ms=2,cam=(0,90),legend=false))

## RK example

#initializing data structs
ts_RK = Time_struct()
@pack! ts_RK = t0,T,dt

arrks_RK = AdjRRK_struct()
@pack! arrks_RK = f,df,η,∇η
@pack! arrks_RK = u0
arrks_RK.u0_lin = u0
arrks_RK.return_time = true
arrks_RK.return_Δη = true

#running forward solver
RK_solver!(arrks_RK,ts_RK,rk4)
Q_RK = arr_to_tuple(arrks_RK.u[:,end])
Q_RK = (x->reshape(x,length(rd.r),md.K)).(Q_RK)
# gr(aspect_ratio=1,legend=false)
display(scatter((x->rd.Vp*x).((x,Q_RK[1])),zcolor=rd.Vp*Q_RK[1],msw=0,ms=2,cam=(0,90),legend=false))
display(plot(ts_RK.t,arrks_RK.Δη,legend=false))

arrks_RK.uT_adj = arrks_RK.u[:,end]

#running linearized solver
RK_solver!(arrks_RK,ts_RK,rk4;lin=true)
W_RK = arr_to_tuple(arrks_RK.u_lin[:,end])
W_RK = (x->reshape(x,length(rd.r),md.K)).(W_RK)
display(scatter((x->rd.Vp*x).((x,W_RK[1])),zcolor=rd.Vp*W_RK[1],msw=0,ms=2,cam=(0,90),legend=false))

#running adjoint solver
RK_solver!(arrks_RK,ts_RK,rk4;adj=true)
Z_RK = arr_to_tuple(arrks_RK.u_adj[:,1])
Z_RK = (x->reshape(x,length(rd.r),md.K)).(Z_RK)
display(scatter((x->rd.Vp*x).((x,Z_RK[3])),zcolor=rd.Vp*Z_RK[3],msw=0,ms=2,cam=(0,90),legend=false))


# #running deriv test on lin RK
# arrks.return_time = false
# arrks.return_Δη = false
# arrks.u0_lin = randn(size(u0))
# Nref = 3
# h0 = 2^(-6)
# arrks_h = AdjRRK_struct(arrks)
# @pack! arrks_h = f,df
# errs,rate,h = AdjRRK.derv_test!(RK_solver!,arrks,arrks_h,ts,rk4,h0,Nref)
# @show rate

## RRK example

#initializing data structs
ts_RRK = Time_struct()
@pack! ts_RRK = t0,T,dt

arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη
@pack! arrks_RRK = u0
arrks_RRK.u0_lin = u0
arrks_RRK.return_time = true
arrks_RRK.return_Δη = true

#running forward solver
RRK_solver!(arrks_RRK,ts_RRK,rk4)
Q_RRK = arr_to_tuple(arrks_RRK.u[:,end])
Q_RRK = (x->reshape(x,length(rd.r),md.K)).(Q_RRK)
display(scatter((x->rd.Vp*x).((x,Q_RRK[1])),zcolor=rd.Vp*Q_RRK[1],msw=0,ms=2,cam=(0,90),legend=false))
display(plot(ts_RRK.t,arrks_RRK.Δη,legend=false))
display(plot(ts_RRK.t,arrks_RRK.γ,legend=false))

arrks_RRK.uT_adj = arrks_RRK.u[:,end]

#running linearized solver
RRK_solver!(arrks_RRK,ts_RRK,rk4;lin=true)
W_RRK = arr_to_tuple(arrks_RRK.u_lin[:,end])
W_RRK = (x->reshape(x,length(rd.r),md.K)).(W_RRK)
display(scatter((x->rd.Vp*x).((x,W_RRK[1])),zcolor=rd.Vp*W_RRK[1],msw=0,ms=2,cam=(0,90),legend=false))

#running adjoint solver
RRK_solver!(arrks_RRK,ts_RRK,rk4;adj=true)
Z_RRK = arr_to_tuple(arrks_RRK.u_adj[:,1])
Z_RRK = (x->reshape(x,length(rd.r),md.K)).(Z_RRK)
display(scatter((x->rd.Vp*x).((x,Z_RRK[1])),zcolor=rd.Vp*Z_RRK[1],msw=0,ms=2,cam=(0,90),legend=false))
