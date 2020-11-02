module AdjRRK

using LinearAlgebra
using UnPack

ENT_TOL = 1e-13
IPT_TOL = 1e-13
DRV_TOL = 5e-2
CNV_TOL = 0.5

struct RKs
    stages::Int64
    b::Array{Float64,1}
    a::Array{Float64,1}
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

    AdjRRK_struct() = new()
end

# coefficients for Heun's method
s_rk2 = 2
b_rk2 = [1/2,1/2]
a_rk2 = [1]
rk2 = RKs(s_rk2,b_rk2,a_rk2)

# coefficients for classic RK4 method
s_rk4 = 4
b_rk4 = [1/6,1/3,1/3,1/6]
a_rk4 = [1/2,1/2,1]
rk4 = RKs(s_rk4,b_rk4,a_rk4)

include("RK_code.jl")
include("RRK_code.jl")
include("test_code.jl")

export RK_solver!, IDT_solver!, RRK_solver
export RKs, rk2, rk4
export AdjRRK_struct, Time_struct, cp_ops
# export derv_test, IP_test

end
