include("KrylovSubspace.jl")

type GMRES{n, T<:Number, TBoundary <: Boundary_Updator} <: GridMatrixIter_one
  A :: Linear_Localized_Functor
  OL_inner :: NTuple{n, Int64}
  pArnoldi_iterator :: ArnoldiIterator{n,T, KrylovSubspace{n,T, Linear_Localized_Functor, TBoundary}}
  function GMRES( A_1 :: Localized_Functor, OL_inner :: NTuple{n,Int64}, ytarget_expand :: Array{T,n}, PBoundary :: TBoundary , param ::Tuple{Vararg{Array{Float64,n}}}, global_param :: Tuple, maximum_order :: Int64)
    OL_bound = PBoundary.OL_bound

    ## A == 0, when inner cube exceeds the boundary
    A = Linear_Localized_Functor(OL_bound, OL_inner, A_1, size(ytarget_expand), param, global_param)
    Subspace = KrylovSubspace(A, ytarget_expand, PBoundary)
    pArnoldi_iterator = ArnoldiIterator(Subspace, maximum_order)
    new(A, OL_inner, pArnoldi_iterator)
  end
end

function solve_expand!(p :: GMRES; iter_times :: Int64 = p.pArnoldi_iterator.maximum_order)
  ## H H^T q = H b
  istep = p.pArnoldi_iterator.istep
  for ii = (istep+1):iter_times
    Forward!(p.pArnoldi_iterator)
  end

  Hessen = p.pArnoldi_iterator.Hessenberg_Upper[1:istep, 1:(istep+1)]
  HHT = Hessen'*Hessen

  ## b = [vecnorm(V_core), 0-000000(n zeros)]
  #= b = fill(0.0, istep+1); =# 
  b_norm = vecnorm(p.Subspace.V_core)
  HTb = fill(0.0, istep)
  Hb[1] = Hessen[1,1] * b_norm;

  projx_s = HHT\Hb
  @assert length(projx_s) == istep

  eigns_expand = p.pArnoldi_iterator.eign_expand[1:istep]
  
  value = fill(0.0 , p.PBoundary.size)
  for ii = 1:istep
    value += eigns_expand[ii] * projx_s[ii]
  end

  error_by_eigns = H*projx_s 
  error_by_eigns[1] -= b_norm

  return (value,error_by_eigns)
end

function Forward!(p :: GMRES)
  Forward!(p.pArnoldi_iterator)
end
