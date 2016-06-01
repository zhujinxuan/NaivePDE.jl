type SOR_AOE{T <: Localized_Functor,n} <: GridMatrixIter_one
  A :: T
  OL_bound :: NTuple{n, Int64}
  OL_inner :: NTuple{n, Int64}
  Range_of_meaningless_y :: NTuple{n, Array{Int64,1}}
end
export SOR_AOE

function get_all_Ainearby{T<: Localized_Functor,n}(
  OL_bound :: NTuple{n,Int64} 
  , OL_inner :: NTuple{n, Int64}
  , Range_of_meaningless_y ::  NTuple{n, Array{Int64,1}}
  , p :: T
  , x :: Array{Float64,n}
  , param :: Tuple{Vararg{Array{Float64,n}}}
  , global_param :: Tuple = ()
  ; iseq :: Array{Int64,1} = Array(Int64,0)
  , delta_x :: Float64 = 0.5)

  irank = map((1:n...)) do ii
    size(x,ii) - 2OL_bound[ii]
  end

  updatepoint_range = ( 0 == length(iseq)) ? collect(1:reduce(*,irank)) : iseq

  inner_rank = map(v->2v+1, OL_inner)
  boundary_rank = map(v->2v+1, OL_bound)
  Ai = fill( fill(NaN, inner_rank), length(updatepoint_range))

  x1background = fill(0.0, size(x))
  x2background = fill(0.0, size(x))

  for (ilogger, ipoint) = enumerate(updatepoint_range)

    center_sub_outer =  begin
      core_sub = ind2sub(irank,ipoint)
      map(v-> core_sub[v] + OL_bound[v], (1:n...))
    end
    x1background[center_sub_outer...] = delta_x

    Ai_set = fill(0.0, inner_rank)

    for ipoint_inner = (1:reduce(*, inner_rank))
      (idx, check_meaningless_y) = begin
        local sub_inner = ind2sub(inner_rank, ipoint_inner)
        local outter_sub = map(v->sub_inner[v]-OL_inner[v]-1+center_sub_outer[v], (1:n...))
        (outter_sub,
        mapreduce( v->in(outter_sub[v], Range_of_meaningless_y[v]), |, (1:n...)))
      end

      if (check_meaningless_y)
        Ai_set[ipoint_inner] = 0.0
      else
        Ax1 = eval(p, x1background, param, idx, global_param)
        Ax2 = eval(p, x2background, param, idx, global_param)
        Ai_set[ipoint_inner] = (Ax1-Ax2)/delta_x
      end
    end

    x1background[center_sub_outer...] = 0.0

    Ai[ilogger] = Ai_set
  end
  return Ai
end


function eval!{T<: Localized_Functor,n}(
  Ja :: SOR_AOE{T,n}
  , x :: Array{Float64,n}
  , y :: Array{Float64,n}
  , param :: Tuple{Vararg{Array{Float64,n}}}
  , global_param :: Tuple = ()
  ; alpha = 1.0, iseq :: Array{Int64,1} = Array(Int64,0)
  , delta_x :: Float64 = 0.5)

  OL_bound = Ja.OL_bound
  OL_inner = Ja.OL_inner
  Range_of_meaningless_y = Ja.Range_of_meaningless_y

  irank = map((1:n...)) do ii
    size(x,ii) - 2OL_bound[ii]
  end

  Ais =  get_all_Ainearby(OL_bound, OL_inner, Range_of_meaningless_y, Ja.A, x, param, global_param; iseq = iseq, delta_x = delta_x)

  updatepoint_range = ( 0 == length(iseq)) ? collect(1:reduce(*,irank)) : iseq

  rand_updatepoint_range = begin
    local rand_benchmark = rand(size(updatepoint_range))
    local rand_ref = sort(collect(1:length(updatepoint_range)), by=v->rand_benchmark[v])
    map(v->(v,updatepoint_range[v]), rand_ref)
  end


  for(ilogger,ipoint)  = rand_updatepoint_range
    Ai = Ais[ilogger]

    center_sub_outer =  begin
      core_sub = ind2sub(irank,ipoint)
      map(v-> core_sub[v] + OL_bound[v], (1:n...))
    end

    outter_subrange = map(v->(center_sub_outer[v]+(-OL_inner[v]:OL_inner[v])), (1:n...))
    ylocal = y[outter_subrange...]

    inner_rank = map(v->2v+1, OL_inner)
    boundary_rank = map(v->2v+1, OL_bound)

    inner_update_range = collect(1:reduce(*,inner_rank))

    Ax = fill(0.0, inner_rank)

    for ipoint_inner = (1:reduce(*, inner_rank))

      (idx, check_meaningless_y) = begin
        local sub_inner = ind2sub(inner_rank, ipoint_inner)
        local outter_sub = map(v->sub_inner[v]-OL_inner[v]-1+center_sub_outer[v], (1:n...))
        (outter_sub,
        mapreduce( v->in(outter_sub[v], Range_of_meaningless_y[v]), |, (1:n...)))
      end

      if(check_meaningless_y)
        Ax[ipoint_inner] = 0.0 
      else
        Ax[ipoint_inner] = eval(Ja.A, x, param, idx, global_param)
      end
    end

    x[center_sub_outer...] += alpha*sum((ylocal - Ax).* Ai)/sum(Ai.* Ai)
    x[center_sub_outer...] += alpha*sum((ylocal - Ax).* Ai)/sum(Ai.* Ai)
  end
  return x
end





