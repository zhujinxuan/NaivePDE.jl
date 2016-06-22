type Linear_Localized_Functor{n, TNumber <: Number} <: Localized_Functor
  OL_bound :: NTuple{n, Int64}
  OL_inner :: NTuple{n, Int64}
  size     :: NTuple{n, Int64}
  meaningless_subs :: NTuple{n, Array{Int64,1}}
  Ai :: Array{Array{TNumber, n}, n }
  A0 :: Array{TNumber,n}
end
export Linear_Localized_Functor

# Do I need PBoundary here?
function Linear_Localized_Functor{n, TNumber <: Number}(
  PBoundary :: Boundary_Updator
  , OL_Local :: NTuple{n, Int64}
  , p :: Localized_Functor
  , param :: Tuple{Vararg{Array{Float64,n}}}
  , global_param :: Tuple, test_value :: TNumber
  ; delta_x :: TNumber = 0.5
  )

  OL_bound = PBoundary.OL_bound
  for ii in 1:n 
    @assert OL_Local[ii] <= OL_bound[ii]
  end


  Local_cube_size = map(v->2v+1, OL_Local)

  Ai = Array(Array{TNumber,n}, ExpandSize(PBoundary))
  A0 = Array(TNumber, ExpandSize(PBoundary))


  x1background = fill(0.0,ExpandSize(PBoundary))
  Core_start = map(v->v+1, OL_bound)
  Core_stop = list_eval(OL_bound, CoreSize(PBoundary))

  for Cartesian_Outter_Point in CartesianRange(CartesianIndex(Core_start), CartesianIndex(Core_stop))
    A0[Cartesian_Outter_Point] = eval(p, x1background, param, idx, global_param)
  end

  for Cartesian_Outter_Point in CartesianRange(CartesianIndex(Core_start), CartesianIndex(Core_stop))
    Ai_set = Array(TNumber,Local_cube_size)

    for Cartesian_Inner_Point in CartesianRange(size(Ai_set))
      Cartesian_idx = Cartesian_Inner_Point + Cartesian_Outter_Point - CartesianIndex(map(v->1+v, OL_Local))
      x1background[Cartesian_idx] = delta_x
      delta_y = eval(p, x1background, param, Cartesian_Outter_Point.I, global_param)
      Ai_set[Cartesian_Inner_Point] = (delta_y-A0[Cartesian_idx])/delta_x 
      x1background[Cartesian_idx] = 0.0
    end
    Ai[Cartesian_Outter_Point] = Ai_set
  end

  Linear_Localized_Functor(OL_bound, OL_inner, outter_size, Range_of_meaningless_y, Ai, A0)
end

function eval{T<: Number, n }(
  p :: Linear_Localized_Functor{n, T}
  , xbackground :: Array{T,n}
  , param :: Tuple{Vararg{Array{T,n}}}
  , idx :: NTuple{n,Int64}
  ,  global_param :: Tuple
  )
  Range_of_meaningless_y = p.meaningless_subs
  local_range = list_eval(idx, p.OL_inner) do v1,c1
    v1+(-c1:c1)
  end
  return dot(collect(xbackground[local_range...]), collect(p.Ai[idx...]))
end
export eval

