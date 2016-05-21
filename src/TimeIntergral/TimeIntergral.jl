## Time Intergral Methods designed for non-parallel programme

abstract TimeIntergral_onestep
abstract ModelParams 
export ModelParams


function TimeForward!( p :: ModelParams, t :: Float64, dt :: Float64)
  error("No function for abstract types")
end


function fTrend( data:: Tuple, p:: ModelParams)
  error("No function for abstract types")
end

function eval!(tf :: TimeIntergral_onestep, p :: ModelParams,
  data :: Tuple, trends :: Tuple,
  t :: Float64, dt :: Float64)
  error("No function for abstract types")
end
export TimeForward

include("AdamsBashforth.jl")
include("RK4.jl")
include("DatasTrends.jl")
