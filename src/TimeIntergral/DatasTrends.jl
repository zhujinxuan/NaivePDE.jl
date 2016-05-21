immutable DatasTrends{norder} <: TimeIntergral_onestep
  n :: Int64
  function DatasTrends()
    @assert eltype(norder) <: Int64
    @assert norder >0
    new(norder)
  end
end
export DatasTrends


function eval!(tf :: DatasTrends{2}, p :: ModelParams,
  data :: Tuple, trends :: Tuple,
  t :: Float64, dt :: Float64)
  ## Trend on y(n)
  ytrend0 = fTrend(data[1],p)
  ## Trend on y(n-1)
  ytrendp = trends
  (y0,yp) = data
  yn = map((1:length(y0)...)) do ii
    correction1 = y0[ii] - yp[ii] - ytrend0[ii]*dt
    correction2 = dt*(ytrend0[ii] - ytrendp[ii])
    y0[ii] + ytrend0[ii]*dt - 5*correction1 + 2*correction2
  end
  return (p,(yn,y0),ytrend0)
end
