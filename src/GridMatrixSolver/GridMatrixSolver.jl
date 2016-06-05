#= abstract GridMatrixIter_one =#
abstract Localized_Functor
export Localized_Functor

#= function eval{T<: Number, n}( =#
#=   A :: Localized_Functor =#
#=   , xbackground :: Array{T,n} =#
#=   , param :: Tuple{Vararg{Array{T,n}}} =#
#=   , indx :: NTuple{n,Int64} =#
#=   ,  global_param :: Tuple =#
#=   ) =#
#=   error("No Function for abstract types") =#

#= end =#

abstract GridMatrixSolver

function solve(p :: GridMatrixSolver, iter_times :: Int64, alpha :: Float64 = 1.0)
  error("No Function for abstract types")
end
export solve
function Forward!(p :: GridMatrixSolver, alpha :: Float64)
  error("No Function for abstract types")
end
export Forward!

include("linear_convertor.jl")
include("Jacobi.jl")
include("SOR.jl")
#= include("SOR_AOE.jl") =#
include("GMRES.jl")


#= include("Error_Checker.jl") =#
