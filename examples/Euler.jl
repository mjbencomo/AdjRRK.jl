# This code runs some numerical experiments and generates .mat files for adjoint
# RRK paper. In particular, we run derivative tests, verifying linearization of
# RRK, on 1D compressible Euler equations.

using Plots
using UnPack
using LinearAlgebra
using SparseArrays
using Random
using LaTeXStrings
using MAT

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

# Flags
write_mat = true  # Want to output .mat file?
make_plot = true # Want to output plots?

if write_mat
    file = matopen("/Users/mariobencomo/Desktop/Research/AdjRRK_paper/figs/RRK_Euler.mat","w")
end

N = 3
K1D = 16
T = 1.5 #0.25
t0 = 0
CN = (N+1)^2/2  # estimated trace constant
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
w = 25*2
rho = @. 1 + .5*exp(-w*(x-.1)^2)
#rho = @. 2.5*(abs.(x)<1/2) + 2*(abs.(x)>=1/2)
#rho = @. 1+ (x<0)
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
Phg = kron(diagm(vec(1 ./ md.J[1,:])),Ph)


QgTr = sparse(transpose(Qg))
global_SBP_ops = (Qg,QgTr,Bg,Vhg,Phg,VhPg,Vqg)

function LF_dissipation(rhoL,rhouL,EL,rhoR,rhouR,ER)
    QL = rhoL,rhouL,EL
    QR = rhoR,rhouR,ER
    λL = abs.(wavespeed_1D(QL...))
    λR = abs.(wavespeed_1D(QR...))
    λavg =.5*(λL+λR)
    #λavg = 0
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

################################################################################
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

Q0 = vec.(Q0)
u0 = tuple_to_arr(Q0)


#==============================================================================#
# Derivative tests

arrks_RRK = AdjRRK_struct()
@pack! arrks_RRK = f,df,η,∇η,Hη
@pack! arrks_RRK = u0
arrks_RRK.return_time = true
arrks_RRK.return_Δη = true

arrks_h = AdjRRK_struct()
@pack! arrks_h = f,df,η,∇η,Hη

ts_RRK = Time_struct()
@pack! ts_RRK = t0,T

Nref = 9
h0 = 2^(-12)

Random.seed!(05052021)
arrks_RRK.u0_lin = randn(size(u0))


# CFL=1 case ------------------------------------------------------------------#
CFL=1
dt = CFL / (CN*K1D)
ts_RRK.dt = dt

# with proper linearization
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = false
errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# dt*-const case
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = true
errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# γ-const case
arrks_RRK.γ_cnst = true
arrks_RRK.dt_cnst = true
errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

if make_plot
    labels = ["proper lin." "Δt*-const" "γ-const"]
    markers = [:circle :square :star5]

    errors = [errs_RRK,errs_RRK_dt0,errs_RRK_γ0]
    plot(title="FD error (CFL=1)",
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
    #     title="FD convergence rate (CFL=1)",
    #     xlabel="refinement index",
    #     ylabel="rate",
    #     label=labels,
    #     marker=markers,
    #     xaxis=:flip)
    # display(plot!())
end

if write_mat
    write(file,"errs_RRK_CFL1",errs_RRK)
    write(file,"errs_RRK_CFL1_dt0",errs_RRK_dt0)
    write(file,"errs_RRK_CFL1_g0",errs_RRK_γ0)
end

# CFL=1.5 case ------------------------------------------------------------------#
CFL=1.5
dt = CFL / (CN*K1D)
ts_RRK.dt = dt

# with proper linearization
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = false
errs_RRK,rate_RRK,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# dt*-const case
arrks_RRK.γ_cnst = false
arrks_RRK.dt_cnst = true
errs_RRK_dt0,rate_RRK_dt0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

# γ-const case
arrks_RRK.γ_cnst = true
arrks_RRK.dt_cnst = true
errs_RRK_γ0,rate_RRK_γ0,h = AdjRRK.derv_test!(RRK_solver!,arrks_RRK,arrks_h,ts_RRK,rk4,h0,Nref)

if make_plot
    labels = ["proper lin." "Δt*-const" "γ-const"]
    markers = [:circle :square :star5]

    errors = [errs_RRK,errs_RRK_dt0,errs_RRK_γ0]
    plot(title="FD error (CFL=1.5)",
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
    #     title="FD convergence rate (CFL=1.5)",
    #     xlabel="refinement index",
    #     ylabel="rate",
    #     label=labels,
    #     marker=markers,
    #     xaxis=:flip)
    # display(plot!())
end

if write_mat
    write(file,"errs_RRK_CFL2",errs_RRK)
    write(file,"errs_RRK_CFL2_dt0",errs_RRK_dt0)
    write(file,"errs_RRK_CFL2_g0",errs_RRK_γ0)
end

if write_mat
    write(file,"h_derv_test",h)
    close(file)
end
