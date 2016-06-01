type Linear_Localized_Functor{n, TNumber <: Number} <: Localized_Functor
  OL_bound :: NTuple{n, Int64}
  OL_inner :: NTuple{n, Int64}
  size     :: NTuple{n, Int64}
  meaningless_subs :: NTuple{n, Array{Int64,1}}
  Ai :: Array{Array{TNumber, n}, n }
  A0 :: Array{TNumber,n}
  function Linear_Localized_Functor(
    OL_bound :: NTuple{n,Int64} 
    , OL_inner :: NTuple{n, Int64}
    , p :: Localized_Functor
    , outter_size :: NTuple{n, Int64}
    , param :: Tuple{Vararg{Array{Float64,n}}}
    , global_param :: Tuple 
    ; delta_x :: Float64 = 0.5
    , Range_of_meaningless_y ::  NTuple{n, Array{Int64,1}} = list_eval(OL_inner, outter_size) do v1,v2
      cat(1,collect(1:v1),collect(v2-v1+(1:v1)))
    end
    )

    isize = map((1:n...)) do ii
      outter_size - 2OL_bound[ii]
    end

    inner_size = map(v->2v+1, OL_inner)
    boundary_size = map(v->2v+1, OL_bound)
    Ai = Array(Array{TNumber,n}, outter_size)
    A0 = Array(TNumber, outter_size)

    x1background = fill(0.0,outter_size)

    for Cartesian_Outter_Point in eachindex(x1background)
      #= iI = Cartesian_Inner_Point =#
      whether_meanless = list_reduce(Cartesian_Outter_Point.I,Range_of_meaningless_y) do x, rr
        in(x,rr)
      end
      if( reduce(|, whether_meanless) )
        A0[Cartesian_Outter_Point] = 0.0
      else
        idx = Cartesian_Outter_Point.I
        A0[Cartesian_Outter_Point] = eval(p, x1background, param, idx, global_param)
      end
    end

    for Cartesian_Outter_Point in CartesianRange(CartesianIndex(map(x->1+x,OL_bound)),CartesianIndex(list_eval(+,isize, OL_bound)))
      Ai_set = Array(TNumber,inner_size)
      x1background[Cartesian_Outter_Point] = delta_x
      for Cartesian_Inner_Point in eachindex(Ai_set)
        Cartesian_idx = Cartesian_Inner_Point + Cartesian_Outter_Point - CartesianIndex(map(v->1+v, OL_inner))
        idx = Cartesian_idx.I

        whether_meanless = list_reduce(idx,Range_of_meaningless_y) do x, rr
          in(x,rr)
        end
        if (reduce(|, whether_meanless))
          Ai_set[Cartesian_Inner_Point] = 0.0
        else
          delta_y = eval(p, x1background, param, idx, global_param)
          Ai_set[Cartesian_Inner_Point] = (delta_y-A0[Cartesian_idx])/delta_x
        end
      end
      x1background[Cartesian_Outter_Point] = 0.0
      Ai[Cartesian_Outter_Point] = Ai_set
    end

    new(OL_bound, OL_inner, outter_size, Range_of_meaningless_y, Ai, A0)
  end
end

function eval{T<: Number, n, }(
  p :: Linear_Localized_Functor{n, T}
  , xbackground :: Array{T,n}
  , param :: Tuple{Vararg{Array{T,n}}}
  , indx :: NTuple{n,Int64}
  ,  global_param :: Tuple
  )
  whether_meanless = list_reduce(idx,Range_of_meaningless_y) do x, rr
    in(x,rr)
  end
  if (reduce(|,whether_meanless))
    return convert(T,0.0)
  else
    local_range = list_eval(idx, OL_inner) do v1,c1
      v1+(-c1:c1)
    end
    return ((xbackground[local_range...] .* p.Ai[idx]) +Ax[idx])
  end
end

