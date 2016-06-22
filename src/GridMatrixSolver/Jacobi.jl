type Jacobi_Solver{n, TNumber <: Number, TBoundary <: Boundary_Updator} <: GridMatrixSolver
  A :: Linear_Localized_Functor{n, TNumber}
  OL_inner :: NTuple{n, Int64}
  PBoundary :: TBoundary
  target :: Array{TNumber, n}
  test_value_expand :: Array{TNumber,n}
  istep :: Int64
end

function Jacobi_Solver{n, TNumber <: Number, TBoundary <: Boundary_Updator}(
  A_1 :: Localized_Functor, OL_inner :: NTuple{n,Int64}
  , param ::Tuple{Vararg{Array{Float64,n}}}, global_param :: Tuple
  , PBoundary :: TBoundary
  , ytarget_expand :: Array{TNumber,n}
  ; test_initial :: Array{TNumber,n} = ytarget_expand
  )
  A = Linear_Localized_Functor(PBoundary, OL_inner, A_1, size(ytarget_expand), param, global_param, test_initial[1])
  Jacobi_Solver(A, OL_inner, PBoundary, ytarget_expand, test_initial, 0) 
end

export Jacobi_Solver


function solve!(p :: Jacobi_Solver, iter_times :: Int64, alpha :: Float64 = 1.0)
  for ii in (p.istep+1):iter_times
    Forward!(p, alpha)
  end

end

function Forward!(p :: Jacobi_Solver, alpha :: Float64)
  p.istep +=1
  istep = p.istep
  PBoundary = p.PBoundary
  OL_bound = PBoundary.OL_bound
  (OL_bound_start, OL_bound_end) =  Boundary_Start_End_subs(PBoundary)
  test_value = p.test_value_expand
  changes_x = zeros(test_value)

  for Cartesian_Inner_Idx in CartesianRange(CartesianIndex(OL_bound_start), CartesianIndex(OL_bound_end))
    idx = Cartesian_Inner_Idx.I
    Ax = eval(p.A, test_value, (), idx, ())
    Ai_c = get_Ai_center(p.A, idx)
    ytarget = p.target[Cartesian_Inner_Idx]
    changes_x[Cartesian_Inner_Idx] = alpha * (ytarget-Ax)/Ai_c 
  end

  p.test_value_expand += changes_x

  Update_Boundaries!(PBoundary, p.test_value_expand)
end
