module AdjRRK

using LinearAlgebra
using UnPack

ENT_TOL = 1e-12
IPT_TOL = 1e-12
DRV_TOL = 0.1
CNV_TOL = 0.5

# struct RKs
#     stages::Int64
#     b::Array{Float64,1}
#     a::Array{Float64,1}
# end

struct RK_struct
    stages::Int64
    b::Array{Float64,1}
    A::Array{Float64,2}
end

mutable struct Time_struct
    t0::Float64
    T::Float64
    dt::Float64
    Nt::Int64
    dt_corr::Float64
    t::Array{Float64,1}
    Time_struct() = new(0,0,0,0,0)
end

mutable struct AdjRRK_struct
    #flags
    return_time::Bool
    return_Δη::Bool
    γ_cnst::Bool
    dt_cnst::Bool

    #functions
    f; df
    η; ∇η; Hη

    #fields
    u::Array{Float64,2}
    u_lin::Array{Float64,2}
    u_adj::Array{Float64,2}

    #initial/final values
    u0::Array{Float64,1}
    u0_lin::Array{Float64,1}
    uT_adj::Array{Float64,1}

    #other
    Δη::Array{Float64,1}
    γ::Array{Float64,1}

    AdjRRK_struct() = new(false,false,false,false)
    AdjRRK_struct(arrks::AdjRRK_struct) = new(arrks.return_time,arrks.return_Δη,arrks.γ_cnst,arrks.dt_cnst)
end

# coefficients for Heun's method
s_rk2 = 2
A_rk2 = [0 0; 1 0]
b_rk2 = [1/2,1/2]
rk2 = RK_struct(s_rk2,b_rk2,A_rk2)

# coefficients for RK3 from Shu & Osher 1988
s_rk3 = 3
A_rk3 = [0 0 0; 1 0 0; 1/4 1/4 0]
b_rk3 = [1/6; 1/6; 2/3]
rk3 = RK_struct(s_rk3,b_rk3,A_rk3)

# coefficients for classic RK4 method
s_rk4 = 4
A_rk4 = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0]
b_rk4 = [1/6,1/3,1/3,1/6]
rk4 = RK_struct(s_rk4,b_rk4,A_rk4)


include("RK_code.jl")
include("RRK_code.jl")
include("test_code.jl")

export RK_solver!, IDT_solver!, RRK_solver!
export RKs, rk2, rk4, rk3, RK_struct
export AdjRRK_struct, Time_struct, cp_ops

end
