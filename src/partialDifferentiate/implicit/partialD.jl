immutable StencilsnN{n_inner, n_outer} <: PartialDImplicitStrategy
  inner :: Int64
  outer :: Int64
  function StencilsnN()
    @assert eltype(n_inner) <: Int64
    @assert eltype(n_outer) <: Int64
    @assert n_inner <= n_outer
    new(n_inner, n_outer)
  end
end
export StencilsnN


# Defaultly set \Delta x = 1.0
# Defaultly set 3/3 Stencils, for simplify LU decompositon
function eval{n1, n2, n3}(
  pf :: PartialD_OP{PartialD{n1}, StencilsnN{3,3}, NaiveBoundary}, 
  data :: Array{Float64,1}; delta_x :: Float64 = 1.0
  )
  datax = initR(pf, data)
  (a,b,c)= initL(pf, data)
  (l,d,u) = naiveLU(a,b,c)
  return BackwardThomas(d,u, ForwardThomas(l,datax)) / ((delta_x)^(pf.pfunctor.n))
end

function eval{n1, n2, n3}(
  pf :: PartialD_OP{PartialD{n1}, StencilsnN{3,3}, NaivePeriodicShared}, 
  data_unshared :: Array{Float64,1}; delta_x :: Float64 = 1.0
  )
  n_shared = pf.bfunctor.n_shared
  data = cat((1,),data_unshared[end-n_shared+(1:n_shared)], data_unshared, data_unshared[1:n_shared] )
  datax = initR(pf, data)
  (a,b,c)= initL(pf, data)
  (l,d,u) = naiveLU(a,b,c)
  res =  BackwardThomas(d,u, ForwardThomas(l,datax)) / ((delta_x)^(pf.pfunctor.n))
  res[1+n_shared:end-n_shared]
end

include("initRL.jl")
include("naiveLU.jl")
