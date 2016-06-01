immutable KrylovSubspace{n, T <: Number, OPT <: Localized_Functor, TBoundary <:Boundary_Updator}
    A :: OPT     
    V_expand :: Array{T,n}
    V_core :: Array{T,n}
    PBoundary :: TBoundary

    function KrylovSubspace(A :: OPT, b_expand :: Array{T,n}, PBoundary :: TBoundary)
      V_expand = b_expand
      V_core = GetCore(PBoundary, b_expand)
      new(A, V_expand, V_core, PBoundary)
    end
end

function Iterator_Expand{n, T<:Number, TKrylovSubspace <: KrylovSubspace}(p :: TKrylovSubspace, b_expand :: Array{T,n})
  PBoundary = p.PBoundary
  res = Array(T,PBoundary.size)
  OL_bound = PBoundary.OL_bound
  OL_bound_start = CartesianIndex{n}(map(v->v+1, OL_bound))
  OL_bound_end = CartesianIndex{n}(list_eval(-, PBoundary.size, OL_bound))
 
  for Cartesian_idx in CartesianRange{n}(OL_bound_start, OL_bound_end)
    idx = Cartesian_idx.I
    param = ()
    global_param = ()
    res[Cartesian_idx] = eval(p.A, b_expand, param, idx, global_param) 
  end
  Update_Boundaries(PBoundary, res)
  return res
end


abstract OrthodoxIterator

type ArnoldiIterator{n, T <: Number, TKrylovSubspace <: KrylovSubspace} <: OrthodoxIterator
  Subspace :: TKrylovSubspace
  eign_expand :: Vector{Array{T,n}}
  eign_core :: Vector{Array{T,n}}
  istep :: Int64
  maximum_order :: Int64
  Hessenberg_Upper :: Array{Float64,2}

  function ArnoldiIterator(Subspace :: TKrylovSubspace, maximum_order :: Int64)

    eign_expand = Vector{Array{T,n}}(Subspace.maximum_order)
    eign_core = Vector{Array{T,n}}(Subspace.maximum_order)
    eign_core[1] = begin
      local unnormalized_core :: Array{T,n}
      unnormalized_core = GetCore(Subspace.PBoundary, Subspace.V_core)
      unnormalized_core / vecnorm(unnormalized_core)
    end
    eign_expand[1] = Expand(Subspace.PBoundary, eign_core[1])

    Hessenberg_Upper = fill(0.0, Subspace.maximum_order, Subspace.maximum_order-1)
    new(Subspace, eign_expand, eign_core, 0, maximum_order, Hessenberg_Upper)
  end
end

function Forward!{n, T<: Number, TKrylovSubspace <: KrylovSubspace}(p :: ArnoldiIterator{n,T,TKrylovSubspace})

  Subspace = p.Subspace

  p.istep += 1
  istep = p.istep

  new_Aq_core = begin
    local new_Ab_expand = Array(T, Subspace.PBoundary.size)
    GetCore(Subspace.PBoundary, Iterator_Expand(Subspace, p.eign_expand[istep]))
  end

  @inbounds @simd for ii = 1:istep
    p.Hessenberg_Upper[ii,istep] =vecdot(new_Aq_core, p.eign_core[ii])
  end

  new_eign_core = new_Aq_core 
  @inbounds for ii = 1:istep
    new_eign_core -= p.eign_core[ii] * p.Hessenberg_Upper[ii,istep]
  end

  p.Hessenberg_Upper[istep+1, istep]  = vecnorm(new_eign_core)
  new_eign_core = new_eign_core./p.Hessenberg_Upper[istep+1, istep]

  p.eign_core[istep+1] = new_eign_core
  p.eign_expand[istep+1] = Expand(PBoundary, new_eign_core)

  ()
end

export Forward!
