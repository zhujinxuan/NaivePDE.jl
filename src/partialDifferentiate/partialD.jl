abstract PartialDN 
immutable PartialD{order} <: PartialDN  
  n :: Int64
  function PartialD()
    @assert eltype(order) <: Int
    return new(order)
  end
end
export PartialD

abstract PartialDStrategy 
abstract PartialDImplicitStrategy <: PartialDStrategy
abstract PartialDExplicitStrategy <: PartialDStrategy

abstract BoundaryStrategy
immutable NaiveBoundary <: BoundaryStrategy end
immutable NaivePeriodicShared <: BoundaryStrategy 
  n_shared :: Int64
  function NaivePeriodicShared(n :: Int64 = 3)
    @assert n > 0
    new(n)
  end
end
export NaiveBoundary
export NaivePeriodicShared

immutable PartialD_OP{T1 <: PartialDN, T2<: PartialDStrategy, T3 <: BoundaryStrategy} <: oneDoperator
  pfunctor :: T1
  sfunctor :: T2
  bfunctor :: T3
  function PartialD_OP()
    new(T1(),T2(),T3())
  end
  function PartialD_OP(x1 :: T1, x2 :: T2, x3:: T3)
    new(x1,x2,x3)
  end
  function PartialD_OP(x1 :: Tuple, x2 :: Tuple, x3:: Tuple)
    new(T1(x1...), T2(x2...), T3(x3...))
  end
end
export PartialD_OP


include("implicit/partialD.jl")

