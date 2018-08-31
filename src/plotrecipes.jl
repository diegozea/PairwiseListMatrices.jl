# ----------------------------- HELPER FUNCTIONS ----------------------------- #

function _getxyz(plm::PairwiseListMatrix{T,D,TV}) where {T,D,TV}
    names = getlabels(plm)
    names, names, Matrix(plm)
end

function _getxyz(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN}
    names = getlabels(nplm)
    names, names, Matrix(nplm.array)
end

# ---------------------------------- MATRIX ---------------------------------- #

@recipe function plot(plm::PairwiseListMatrix{T,D,TV}) where {T,D,TV}
    seriestype  :=  :heatmap
    yflip --> true
    ratio --> :equal
    _getxyz(plm)
end

@recipe function plot(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN}
    seriestype  :=  :heatmap
    yflip --> true
    ratio --> :equal
    _getxyz(nplm)
end
