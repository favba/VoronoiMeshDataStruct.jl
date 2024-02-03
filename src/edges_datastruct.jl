struct EdgeConnectivity{TI<:Integer}
    """Indices of the vertices that form the edge"""
    vertices::Vector{NTuple{2,TI}}
    """Indices of the cells divided by the edge"""
    cells::Vector{NTuple{2,TI}}
end

struct Edges{TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Edges connectivity data struct"""
    indices::EdgeConnectivity{TI}
    """Edge's x,y,z coordinates"""
    position::VAPos
end