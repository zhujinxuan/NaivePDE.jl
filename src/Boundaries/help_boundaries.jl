abstract Boundary_Updator

function Update_Boundaries!{n, TFloat <: Number}(
  pbound :: Boundary_Updator, data :: Array{TFloat, n})
  error("No Methods for an undefined type")
end

function GetCore{n, TFloat <: Number}(
  pbound :: Boundary_Updator, data :: Array{TFloat, n})
  error("No Function for abstract types")
end

function Expand{n, TFloat <: Number}(
  pbound :: Boundary_Updator, data :: Array{TFloat, n})
  error("No Function for abstract types")
end

function Boundary_Start_End_subs(pbound :: Boundary_Updator)
  OL_bound = pbound.OL_bound
  bound_starts = map(v->1+v, OL_bound)
  bound_ends   = list_eval(-, pbound.size, OL_bound)
  return (bound_starts, bound_ends)
end

include("Periodic.jl")
#= include("Complex.jl") =#

## Some Simple Functions

#= function periodic_makes_boundary!{T1 <: Number,n}( =#
#=   x:: Array{T1,n}, OL :: NTuple{n, Int64}) =#

#=   irank = map((1:n...)) do ii =#
#=     size(x,ii) - 2OL[ii] =#
#=   end =#

#=   xcore = Array(T1,irank) =#
#=   xcore[:] = x[ map( ii-> OL[ii] + (1:irank[ii]) , (1:n...))...] =#

#=   for ipoint = 1:length(x) =#
#=     irr = ind2sub(size(x),ipoint) =#

#=     back_to_core_sub = map((1:n...)) do ii =#
#=       mod(irr[ii] - OL[ii]-1, irank[ii])+1 =#
#=     end =#
#=     x[irr...] = xcore[ back_to_core_sub...] =#
#=   end =#

#=   return (x, xcore) =#

#= end =#
#= export periodic_makes_boundary! =#
  

#= function periodic_expand_boundary{T1 <: Number,n}( =#
#=   xcore:: Array{T1,n}, OL :: NTuple{n, Int64}) =#

#=   expandrank = map((1:n...)) do ii =#
#=     size(xcore,ii) + 2OL[ii] =#
#=   end =#

#=   x = Array(T1, expandrank) =#

#=   for ipoint = 1:length(x) =#
#=     irr = ind2sub(size(x),ipoint) =#

#=     back_to_core_sub = map((1:n...)) do ii =#
#=       mod(irr[ii] - OL[ii]-1, size(xcore,ii))+1 =#
#=     end =#
#=     x[irr...] = xcore[ back_to_core_sub...] =#
#=   end =#

#=   return (x, xcore) =#

#= end =#
#= export periodic_expand_boundary =#


#= function get_core{T1 <: Number,n}( =#
#=   x:: Array{T1,n}, OL :: NTuple{n, Int64}) =#
#=   core_sub_range = map((1:n...)) do ii =#
#=     (OL[ii]+1):(size(x,ii)-OL[ii]) =#
#=   end =#
#=   x1 = x[core_sub_range...] =#
#=   return x1 =#
#= end =#
#= export get_core =#

#= function NaN_expand{T1 <: Number,n}( =#
#=   xcore:: Array{T1,n}, OL :: NTuple{n, Int64}) =#

#=   core_sub_range = map((1:n...)) do ii =#
#=     (OL[ii]+1):(size(x,ii)-OL) =#
#=   end =#

#=   expandrank = map((1:n...)) do ii =#
#=     size(x,ii) + 2OL[ii] =#
#=   end =#

#=   x = Array(T1,expandrank) =#
#=   x[:] = NaN =#
#=   x[core_sub_range...] = xcore =#
  
#= end =#

