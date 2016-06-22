abstract Boundary_Updator{n}

function Update_Boundaries!{n, TFloat <: Number}(
  pbound :: Boundary_Updator{n}, data :: Array{TFloat, n})
  error("No Methods for an undefined type")
end

CoreSize{n}(pbound :: Boundary_Updator{n}) = list_eval(-,pbound.size, map(x->2x, pbound.OL_bound))
ExpandSize{n}(pbound :: Boundary_Updator{n}) = pbound.size

function GetCore{n, TFloat <: Number}(
  pbound :: Boundary_Updator{n}, data :: Array{TFloat, n})
  @assert size(data) == ExpandSize(pbound)
  OL_bound = pbound.OL_bound
  core_subs = map( (1:ndims(data)...)) do ii
    bound_length = OL_bound[ii]
    core_length = size(data,ii) - 2bound_length
    (bound_length + (1:core_length))
  end
  return data[core_subs...]
end

function Expand{n, TFloat <: Number}(
  pbound :: Boundary_Updator{n}, data :: Array{TFloat, n};
  initial_by_ :: Function = zeros)
  @assert size(data) == CoreSize(pbound)

  OL_bound = pbound.OL_bound
  data_expand = initial_by_(TFloat, pbound.size)

  left_register = list_eval(size(data),OL_bound) do core_length, bound_length
    (bound_length+1):(core_length+bound_length)
  end
  data_expand[left_register...] = data
  Update_Boundaries!(pbound, data_expand)
  return data_expand
end

function CoreSubs{n}(pbound :: Boundary_Updator{n})
  OL_bound = pbound.OL_bound
  return list_eval(pbound.size, OL_bound) do xwhole, xbound
    (1+xbound):(xwhole-xbound)
  end
end

function Boundary_Start_End_subs{n}(pbound :: Boundary_Updator{n})
  OL_bound = pbound.OL_bound
  bound_starts = map(v->1+v, OL_bound)
  bound_ends   = list_eval(-, pbound.size, OL_bound)
  return (bound_starts, bound_ends)
end
export Expand, GetCore, Boundary_Start_End_subs

include("Periodic.jl")
include("numeric.jl")
