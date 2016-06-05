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
