# ----------------------------- HELPER FUNCTIONS ----------------------------- #

function _getxyz{T,D,TV}(plm::PairwiseListMatrix{T,D,TV})
    names = getlabels(plm)
    names, names, full(plm)
end

function _getxyz{T,D,TV,DN}(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN})
    names = getlabels(nplm)
    names, names, full(nplm.array)
end

# ---------------------------------- MATRIX ---------------------------------- #

@recipe function plot{T,D,TV}(plm::PairwiseListMatrix{T,D,TV})
    seriestype  :=  :heatmap
    yflip --> true
    ratio --> :equal
    _getxyz(plm)
end

@recipe function plot{T,D,TV,DN}(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN})
    seriestype  :=  :heatmap
    yflip --> true
    ratio --> :equal
    _getxyz(nplm)
end
