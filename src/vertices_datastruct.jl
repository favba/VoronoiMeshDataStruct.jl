struct VertexConnectivity{TI<:Integer}
    """Indices of the edges that meet at the vertex"""
    edges::Vector{NTuple{3,TI}}
    """Indices of the cells that meet at the vertex"""
    cells::Vector{NTuple{3,TI}}
end

integer_precision(::Type{<:VertexConnectivity{T}}) where T = T

Base.getproperty(vertex::VertexConnectivity,s::Symbol) = _getproperty(vertex,Val(s))
_getproperty(vertex::VertexConnectivity,::Val{s}) where s = getfield(vertex,s)
_getproperty(vertex::VertexConnectivity,::Val{:edgesOnVertex}) = getfield(vertex,:edges)
_getproperty(vertex::VertexConnectivity,::Val{:cellsOnVertex}) = getfield(vertex,:cells)

struct VertexBase{TI<:Integer, VAPos<:VecArray{<:Any,1}}
    n::Int
    """Vertices connectivity data struct"""
    indices::VertexConnectivity{TI}
    """Vertex's x,y,z coordinates"""
    position::VAPos
end

integer_precision(::Type{<:VertexBase{T}}) where T = T
float_precision(::Type{<:VertexBase{T,V}}) where {T,V} = TensorsLite._my_eltype(eltype(V))

Base.getproperty(vertex::VertexBase,s::Symbol) = _getproperty(vertex,Val(s))
_getproperty(vertex::VertexBase,::Val{s}) where s = getfield(vertex,s)
_getproperty(vertex::VertexBase,::Val{:nVertices}) where s = getfield(vertex,:n)
_getproperty(vertex::VertexBase,::Val{:edgesOnVertex}) = getfield(vertex,:indices).edges
_getproperty(vertex::VertexBase,::Val{:cellsOnVertex}) = getfield(vertex,:indices).cells
_getproperty(vertex::VertexBase,::Val{:xVertex}) = getfield(vertex,:position).x
_getproperty(vertex::VertexBase,::Val{:yVertex}) = getfield(vertex,:position).y
_getproperty(vertex::VertexBase,::Val{:zVertex}) = getfield(vertex,:position).z

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

Base.getproperty(vertex::VertexInfo,s::Symbol) = _getproperty(vertex,Val(s))

_getproperty(vertex::VertexInfo,::Val{s}) where s = getfield(vertex,s)
_getproperty(vertex::VertexInfo,::Val{:n}) = getfield(vertex,:base).n
_getproperty(vertex::VertexInfo,::Val{:indices}) = getfield(vertex,:base).indices
_getproperty(vertex::VertexInfo,::Val{:position}) = getfield(vertex,:base).position

for s in (:nVertices,:cellsOnVertex,:edgesOnVertex,:xVertex,:yVertex,:zVertex)
    @eval _getproperty(vertex::VertexInfo,::Val{$(QuoteNode(s))}) = getproperty(getfield(vertex,:base),$(QuoteNode(s)))
end

for (s,nc) in pairs((longitude=:lonVertex, latitude=:latVertex,
                     indexToID=:indexToVertexID, area=:areaTriangle,
                     bdyMask=:bdyMaskVertex, kiteAreas=:kiteAreasOnVertex))

    @eval _getproperty(vertex::VertexInfo,::Val{$(QuoteNode(nc))}) = getfield(vertex,$(QuoteNode(s)))
end