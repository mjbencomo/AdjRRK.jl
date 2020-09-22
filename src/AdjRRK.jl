module AdjRRK

using LinearAlgebra
#using Interpolations

struct RKs
    stages::Int64
    b::Array
    a::Array
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

export RK_solver, RRK_solver
export RKs, rk2, rk4

end
