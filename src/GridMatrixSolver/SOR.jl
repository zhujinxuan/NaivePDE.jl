type SOR_Solver{n, TNumber <: Number, TBoundary <: Boundary_Updator} <: GridMatrixSolver
  A :: Linear_Localized_Functor{n, TNumber}
  OL_inner :: NTuple{n, Int64}
  PBoundary :: TBoundary
  target :: Array{TNumber, n}
  test_value_expand :: Array{TNumber,n}
  istep :: Int64
end

function SOR_Solver{n, TNumber <: Number, TBoundary <: Boundary_Updator}(
  A_1 :: Localized_Functor, OL_inner :: NTuple{n,Int64}
  , param ::Tuple{Vararg{Array{Float64,n}}}, global_param :: Tuple
  , PBoundary :: TBoundary
  , ytarget_expand :: Array{TNumber,n}
  ; test_initial :: Array{TNumber,n} = ytarget_expand
  )
  A = Linear_Localized_Functor(PBoundary, OL_inner, A_1, size(ytarget_expand), param, global_param, test_initial[1])
  SOR_Solver(A, OL_inner, PBoundary, ytarget_expand, test_initial, 0) 
end

export SOR_Solver


function solve!(p :: SOR_Solver, iter_times :: Int64, alpha :: Float64 = 1.0)
  for ii in (p.istep+1):iter_times
    Forward!(p, alpha)
  end
end

function Forward!(p :: SOR_Solver, alpha :: Float64)
  p.istep +=1
  istep = p.istep
  PBoundary = p.PBoundary
  OL_bound = PBoundary.OL_bound
  (OL_bound_start, OL_bound_end) =  Boundary_Start_End_subs(PBoundary)
  test_value = p.test_value_expand
  for Cartesianrr_sliced in slice_an_range(CartesianRange(CartesianIndex(OL_bound_start), CartesianIndex(OL_bound_end)))
    for Cartesian_Inner_Idx in Cartesianrr_sliced
      idx = Cartesian_Inner_Idx.I
      Ax = eval(p.A, test_value, (), idx, ())
      Ai_c = get_Ai_center(p.A, idx)
      ytarget = p.target[Cartesian_Inner_Idx]
      p.test_value_expand[Cartesian_Inner_Idx] += alpha * (ytarget-Ax)/Ai_c 
    end
    Update_Boundaries!(PBoundary, p.test_value_expand)
  end
end

export SOR_Solver

function slice_an_range{n}( OL_range :: CartesianRange{CartesianIndex{n}})
  bound_start = OL_range.start.I
  bound_stop = OL_range.stop.I
  mid_subs = list_eval(bound_start, bound_stop) do x1,x2
    round(Int64, (x1 + x2)/2)
  end
  
  res = Array{CartesianRange{CartesianIndex{n}}}(map(v->2, (1:n...)))

  for Carteidx in CartesianRange(size(res))
    res_start= map((1:n...)) do ii
      (1==Carteidx[ii] ) ? bound_start[ii] : (mid_subs[ii]+1)
    end
    res_stop = map((1:n...)) do ii
      (1==Carteidx[ii] ) ? mid_subs[ii] : (bound_stop[ii])
    end
    
    res[Carteidx] = CartesianRange(CartesianIndex(res_start), CartesianIndex(res_stop))
  end
  return collect(res)
end
