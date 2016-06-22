function dot_boundary{n, TNumber<: Number}(pb :: Boundary_Updator{n}, a1 :: Array{TNumber,n}, a2 :: Array{TNumber,n})
  return dot_core(pb, GetCore(pb, a1), GetCore(pb,a2))
end

function dot_core{n, TNumber <: Number}(pb :: Boundary_Updator{n}, a1:: Array{TNumber,n}, a2 :: Array{TNumber,n})
  @fastmath vecdot(a1, a2)
end

function normalize_core!{n, TNumber <: Number}(pb :: Boundary_Updator{n}, a1 :: Array{TNumber,n})
   norm_as = sqrt(dot_core(pb, a1, a1))
   a1[:] /= norm_as
   return (norm_as)
end

function normalize_boundary!{n, TNumber<: Number}(pb :: Boundary_Updator{n}, a1 :: Array{TNumber,n})
  norm_as = sqrt(dot_boundary(pb, a1,a1))
  a1[CoreSubs(a1)...] /= norm_as
  Update_Boundaries!(pb, a1)
  return norm_as
end

function proj_core!{n, TNumber<: Number}(pb :: Boundary_Updator{n}, y :: Array{TNumber,n}, eigns :: Vector{Array{TNumber,n}}, eign_nums :: Int64)
  @inbounds @simd for ii = 1:eign_nums
    @fastmath y[:] -= dot_core(pb, y, eigns[ii])* eigns[ii]
  end
  y
end

