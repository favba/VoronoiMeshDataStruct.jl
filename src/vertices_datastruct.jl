struct VertexConnectivity{TI<:Integer}
    """Indices of the edges that meet at the vertex"""
    edges::Vector{NTuple{3,TI}}
    """Indices of the cells that meet at the vertex"""
    cells::Vector{NTuple{3,TI}}
end

struct Vertices{TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Vertices connectivity data struct"""
    indices::VertexConnectivity{TI}
    """Vertex's x,y,z coordinates"""
    position::VAPos
end