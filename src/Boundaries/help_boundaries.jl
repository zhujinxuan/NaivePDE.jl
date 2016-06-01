abstract Boundary_Updator

type Periodic_Boundaries{n}
  OL_bound :: NTuple{n, Int64}
  size :: NTuple{n, Int64}
end
export Periodic_Boundaries



function Update_Boundaries!{n, TFloat <: Number}(
  pbound :: Periodic_Boundaries{n}, data :: Array{TFloat, n})
  @assert size(data) == pbound.size
  OL_bound = pbound.OL_bound
  full_area_ranges = map(v->collect(1:v),size(data))
  converting_core_full_area_ranges = map( (1:ndims(data)...)) do ii
    whole_length = size(data,ii) 
    bound_length = OL_bound[ii]
    core_length = whole_length - 2bound_length
    cat(1,collect(core_length + (1:bound_length)),
          collect(bound_length + (1:core_length)),
          collect(bound_length+(1:bound_length)))
  end

  left_register = Vector{Vector{Int64}}(ndims(data))
  left_register[:] = collect(full_area_ranges)
  right_register = Vector{Vector{Int64}}(ndims(data))
  right_register[:] = collect(converting_core_full_area_ranges)

  for ii = 1:ndims(data)
    if (OL_bound[ii] > 0 )
      left_register[:]  = collect(full_area_ranges)
      right_register[:] = collect(converting_core_full_area_ranges)

      left_pending1 = 0
      left_pending2 = size(data,ii)-OL_bound[ii]
      left_register[ii] = cat(1,collect(1:OL_bound[ii]),
                                 collect(left_pending2+(1:OL_bound[ii])))

      right_pending1 = size(data,ii)-2OL_bound[ii]
      right_pending2 = OL_bound[ii]
      right_register[ii] = cat(1,collect(right_pending1+(1:OL_bound[ii])),
                               collect(right_pending2+(1:OL_bound[ii])))
      data[left_register...] = data[right_register...]
    end
  end
  return data
end

export Update_Boundaries!


function GetCore{n, TFloat <: Number}(
  pbound :: Periodic_Boundaries{n}, data :: Array{TFloat, n})
  @assert size(data) == pbound.size
  core_subs = map( (1:ndims(data)...)) do ii
    bound_length = OL_bound[ii]
    core_length = size(data,ii) - 2bound_length
    (bound_length + (1:core_length))
  end
  return data[core_subs...]
end

function Expand{n, TFloat <: Number}(
  pbound :: Periodic_Boundaries{n}, data :: Array{TFloat, n})
  @assert size(data) == list_eval(-,pbound.size, map(x->2x, pbound.OL_bound))
  OL_bound = pbound.OL_bound
  right_register = list_eval(size(data),OL_bound) do whole_length, bound_length
    cat(1,collect(whole_length-bound_length+(1:bound_length)),
          collect(1:whole_length),
          collect(1:bound_length))
  end
  data_expand = data[right_register...]
  return data_expand
end

export Expand

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

