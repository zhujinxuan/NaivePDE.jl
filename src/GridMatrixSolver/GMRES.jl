include("KrylovSubspace.jl")

type GMRES{n, TNumber<:Number, TArnoldiIterator<:ArnoldiIterator} <: GridMatrixSolver
  A :: Linear_Localized_Functor{n, TNumber}
  OL_inner :: NTuple{n, Int64}
  pArnoldi_iterator :: TArnoldiIterator
  test_value_expand :: Array{TNumber,n}
end

function GMRES{n, TNumber<:Number, TBoundary <: Boundary_Updator}(
  A_1 :: Localized_Functor, OL_inner :: NTuple{n,Int64}
  , param ::Tuple{Vararg{Array{Float64,n}}}, global_param :: Tuple
  , PBoundary :: TBoundary
  , ytarget_expand :: Array{TNumber,n}
  , maximum_order :: Int64)

  OL_bound = PBoundary.OL_bound
  A = Linear_Localized_Functor(PBoundary, OL_inner, A_1, size(ytarget_expand), param, global_param, ytarget_expand[1])
  Subspace = KrylovSubspace(A, ytarget_expand, PBoundary)
  pArnoldi_iterator = ArnoldiIterator(Subspace, maximum_order)
  GMRES(A, OL_inner, pArnoldi_iterator, ytarget_expand)
end

export GMRES

function solve!(pGMRES :: GMRES, iter_times :: Int64 = pGMRES.pArnoldi_iterator.maximum_order -1)
  # H^T H q = H^T b
  istep = pGMRES.pArnoldi_iterator.istep
  if (istep < iter_times)
    for ii = (istep+1):iter_times
      (0==mod(ii,10)) ? println(ii) : nothing
      Forward!(pGMRES.pArnoldi_iterator)
    end
  end
  istep = pGMRES.pArnoldi_iterator.istep

  Hessen = pGMRES.pArnoldi_iterator.Hessenberg_Upper[1:(istep+1), 1:istep]
  HHT = Hessen'*Hessen

  b_norm = vecnorm(pGMRES.pArnoldi_iterator.Subspace.V_core)
  HTb = collect(Hessen[1,:]) * b_norm

  projx_s = HHT\HTb

  eigns_expand = pGMRES.pArnoldi_iterator.eign_expand[1:istep]
  
  value = fill(0.0 , pGMRES.pArnoldi_iterator.Subspace.PBoundary.size)
  for ii = 1:istep
    value += eigns_expand[ii] * projx_s[ii]
  end

  error_by_eigns = Hessen*projx_s
  error_by_eigns[1] -= b_norm
  pGMRES.test_value_expand = value

  return (value,error_by_eigns)
end

function Forward!(p :: GMRES)
  @assert p.pArnoldi_iterator.istep < p.pArnoldi_iterator.maximum_order
  Forward!(p.pArnoldi_iterator)
end
export Forward!
