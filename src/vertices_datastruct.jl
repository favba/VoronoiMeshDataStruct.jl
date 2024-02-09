struct VertexConnectivity{TI<:Integer}
    """Indices of the edges that meet at the vertex"""
    edges::Vector{NTuple{3,TI}}
    """Indices of the cells that meet at the vertex"""
    cells::Vector{NTuple{3,TI}}
end

integer_precision(::Type{<:VertexConnectivity{T}}) where T = T

struct VertexBase{TI<:Integer, VAPos<:VecArray{<:Any,1}}
    """Vertices connectivity data struct"""
    indices::VertexConnectivity{TI}
    """Vertex's x,y,z coordinates"""
    position::VAPos
end

integer_precision(::Type{<:VertexBase{T}}) where T = T
float_precision(::Type{<:VertexBase{T,V}}) where {T,V} = TensorsLite._my_eltype(eltype(V))

mutable struct VertexInfo{TVertex<:VertexBase,TI,TF}
    const base::TVertex
    longitude::Vector{TF}
    latitude::Vector{TF}
    """Mapping from local array index to global vertex ID"""
    indexToID::Vector{TI}
    """Spherical area of Delaunay triangle"""
    area::Vector{TF}
    """Intersection areas between primal (Voronoi) and dual (triangular) mesh cells"""
    kiteAreas::Vector{NTuple{3,TF}}
    """Indicator of whether a vertex is an interior vertex, a relaxation-zone vertex, or a specified-zone vertex"""
    bdyMask::Vector{TI}

    function VertexInfo(vertex::TVertex) where TVertex<:VertexBase
        return new{TVertex,integer_precision(TVertex),float_precision(TVertex)}(vertex)
    end
end

integer_precision(::Type{<:VertexInfo{TV,TI}}) where {TV,TI} = TI
float_precision(::Type{<:VertexInfo{TV,TI,TF}}) where {TV,TI,TF} = TF

@inline Base.getproperty(vertex::VertexInfo,s::Symbol) = _getproperty(vertex,Val(s))

@inline _getproperty(vertex::VertexInfo,::Val{s}) where s = getfield(vertex,s)
@inline _getproperty(vertex::VertexInfo,::Val{:indices}) = getfield(vertex,:base).indices
@inline _getproperty(vertex::VertexInfo,::Val{:position}) = getfield(vertex,:base).position