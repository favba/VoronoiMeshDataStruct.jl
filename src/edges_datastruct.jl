struct EdgeConnectivity{TI<:Integer}
    """Indices of the vertices that form the edge"""
    vertices::Vector{NTuple{2,TI}}
    """Indices of the cells divided by the edge"""
    cells::Vector{NTuple{2,TI}}
end

integer_precision(::Type{<:EdgeConnectivity{TI}}) where TI = TI

struct EdgeBase{TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Edges connectivity data struct"""
    indices::EdgeConnectivity{TI}
    """Edge's x,y,z coordinates"""
    position::VAPos
end

integer_precision(::Type{<:EdgeBase{TI}}) where TI = TI
float_precision(::Type{<:EdgeBase{T,V}}) where {T,V} = TensorsLite._my_eltype(eltype(V))

struct EdgeVelocityReconstruction{MAX_N_EDGES,TI<:Integer,TF}
    """Number of edges involved in reconstruction of tangential velocity for an edge"""
    nEdges::Vector{TI}
    """IDs of edges involved in reconstruction of tangential velocity for an edge"""
    indices::Vector{VariableLengthIndices{MAX_N_EDGES,TI}}
    """Weights used in reconstruction of tangential velocity for an edge"""
    weights::Matrix{TF}
end

max_n_edges_vel_reconstruction(::Type{<:EdgeVelocityReconstruction{N}}) where N = N

mutable struct EdgeInfo{TEdgeBase,TI<:Integer,TF<:Real,TEdgeVelRecon<:EdgeVelocityReconstruction}
    const base::TEdgeBase
    const velRecon::TEdgeVelRecon
    longitude::Vector{TF}
    latitude::Vector{TF}
    """Mapping from local array index to global edge ID"""
    indexToID::Vector{TI}
    """Spherical distance between vertex endpoints of an edge"""
    dv::Vector{TF}
    """Spherical distance between cells separated by an edge"""
    dc::Vector{TF}
    """Angle between local north and the positive tangential direction of an edge"""
    angle::Vector{TF}
    """Indicator of whether an edge is an interior edge, a relation-zone edge, or a specified-zone edge"""
    bdyMask::Vector{TI}
    """Cartesian components of the vector normal to an edge and tangential to the surface of the sphere"""
    normalVectors::Vec3DArray{TF,1}
    """Weights for cell-centered second derivative, normal to edge, for transport scheme"""
    derivTwo::Array{TF,3}

    function EdgeInfo(edge::TEdgeBase,velRecon::TEdgeVelRecon) where {TEdgeBase<:EdgeBase, TEdgeVelRecon<:EdgeVelocityReconstruction}
        return new{TEdgeBase,integer_precision(TEdgeBase),float_precision(TEdgeBase),TEdgeVelRecon}(edge,velRecon)
    end

end

integer_precision(::Type{<:EdgeInfo{TB,TI}}) where {TB,TI} = TI
float_precision(::Type{<:EdgeInfo{TB,TI,TF}}) where {TB,TI,TF} = TF
max_n_edges_vel_reconstruction(::Type{<:EdgeInfo{TB,TI,TF,TEV}}) where {TB,TI,TF,TEV} = max_n_edges_vel_reconstruction(TEV)