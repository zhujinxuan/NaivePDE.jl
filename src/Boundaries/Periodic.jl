immutable Periodic_Boundaries{n} <: Boundary_Updator{n}
  OL_bound :: NTuple{n, Int64}
  size :: NTuple{n, Int64}
  Periodic_Dims :: Tuple{Vararg{Int64}}
  #= function Periodic_Boundaries{n}(OL_bound :: NTuple{n, Int64}, size :: NTuple{n, Int64},  Periodic_Dims :: Tuple{Vararg{Int64}} =# 
    #= new(OL_bound, size, Periodic_Dims) =#
  #= end =#
end
Periodic_Boundaries{n}(OL_bound :: NTuple{n,Int64}, size :: NTuple{n,Int64}, Periodic_Dims :: Tuple{Vararg{Int64}} = (1:n...)) = Periodic_Boundaries{n}(OL_bound, size, Periodic_Dims)
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

  for ii = pbound.Periodic_Dims
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



