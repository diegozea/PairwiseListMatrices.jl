function Plots._apply_recipe(d::Dict, plm::PairwiseListMatrix; kw...)

    get!(d, :yflip, true)
    get!(d, :linetype, :heatmap)

    # args
    (plm, )
end

function arcdiagram(plm::PairwiseListMatrix)

  n = plm.nelements
  x = 1:n

  plt = scatter(x, zeros(n), m=(10,:orange), leg=false)

  zmin,zmax = extrema(plm)
  grad = ColorGradient(:bluesreds)
  curvecolor(v) = getColorZ(grad, (v-zmin)/(zmax-zmin))

  for i in 1:(n-1)
    for j in i:n
                                                       # i-j
      curve = BezierCurve(P2[(x[i],0), ((x[i]+x[j])/2, mat[i,j]), (x[j], 0)])
      plot!(curve_points(curve), line = (curvecolor(mat[i,j]), 0.5, 2))
    end
  end

  plt
end
