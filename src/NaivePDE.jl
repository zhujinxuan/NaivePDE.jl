module NaivePDE

# package code goes here
include("oneDopeartor.jl")
include("TimeIntergral/TimeIntergral.jl")
include("Boundaries/help_boundaries.jl")
include("GridMatrixSolver/GridMatrixSolver.jl")
include("hypobolicPDE/schemes.jl")

export eval

function list_eval{n}(f :: Function, data :: NTuple{n}...)
  map((1:n...)) do ii
    f( map(vv->(data[vv])[ii], (1:length(data)...))...)
  end
end

function list_reduce{n}(f :: Function, data :: NTuple{n}...)
  map((1:n...)) do ii
    reduce(f, map(vv->(data[vv])[ii], (1:length(data)...)))
  end
end
end # module
