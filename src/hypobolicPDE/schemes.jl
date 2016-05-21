abstract hypobolicSolver_oneD <: oneDoperator

type upwind{n} <: hypobolicSolver_oneD
  f :: Function
end

type LaxWendroff <: hypobolicSolver_oneD
  f :: Function
end

type DiscreteData_oneD 
  value :: Array{Float64,1}
  dx :: Float64
end

function eval!(p :: hypobolicSolver_oneD, data :: DiscreteData_oneD, t :: Float64, dt :: Float64)
  error("No Function for abstract types")
end

function eval(p :: upwind{1}, data :: DiscreteData_oneD , t :: Float64, dt :: Float64)
  y1 = data.value
  dx = data.dx
  fy1 = p.f(y1)
  fy2 = cat(1,fy1[end], fy1[1:end-1] )
  y2 = y1 - (fy1 - fy2) * dt/dx
  return DiscreteData_oneD(y2,dx)
end

function eval(p :: upwind{2}, data :: DiscreteData_oneD , t :: Float64, dt :: Float64)
  y1 = data.value
  dx = data.dx
  fy1 = p.f(y1)
  fy2 = cat(1, fy1[end], fy1[1:end-1])
  fy3 = cat(1, fy1[end-1:end], fy1[1:end-2])
  y2 = y1 - (+fy3 -4fy2 + 3fy1) * dt / (2dx)
  return DiscreteData_oneD(y2,dx)
end

function eval(p :: LaxWendroff, data :: DiscreteData_oneD , t :: Float64, dt :: Float64)
  y1 = data.value
  dx = data.dx
  fy1 = p.f(y1)
  
  yright = (cat(1, fy1[2:end], fy1[1])+ fy1)/2
  yleft = (cat(1, fy1[end], fy1[1:end-1])+ fy1)/2
  yright = yright - (cat(1, fy1[2:end], fy1[1]) - fy1) * dt / (2dx)
  yleft = yleft  - (fy1 - cat(1, fy1[end], fy1[1:end-1])) * dt / (2dx)
  
  y2 = y1 - (p.f(yright) - p.f(yleft)) * dt/dx
  return DiscreteData_oneD(y2,dx)
end
