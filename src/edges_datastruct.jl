struct EdgeConnectivity{TI<:Integer}
    """Indices of the vertices that form the edge"""
    vertices::Vector{NTuple{2,TI}}
    """Indices of the cells divided by the edge"""
    cells::Vector{NTuple{2,TI}}
end

integer_precision(::Type{<:EdgeConnectivity{TI}}) where TI = TI

Base.getproperty(edge::EdgeConnectivity,s::Symbol) = _getproperty(edge,Val(s))
_getproperty(edge::EdgeConnectivity,::Val{s}) where s = getfield(edge,s)
_getproperty(edge::EdgeConnectivity,::Val{:verticesOnEdge}) = getfield(edge,:vertices)
_getproperty(edge::EdgeConnectivity,::Val{:cellsOnEdge}) = getfield(edge,:cells)


struct EdgeBase{S,TI<:Integer, TF<:Real,Tz<:Number}
    n::Int
    """Edges connectivity data struct"""
    indices::EdgeConnectivity{TI}
    """Edge's x,y,z coordinates"""
    position::TensorsLite.VecMaybe2DxyArray{TF,Tz,1}
    onSphere::Val{S}
end

on_a_sphere(::Type{<:EdgeBase{B}}) where B = B
integer_precision(::Type{<:EdgeBase{B,T}}) where {B,T} = T
float_precision(::Type{<:EdgeBase{B,T,TF}}) where {B,T,TF} = TF

Base.getproperty(edge::EdgeBase,s::Symbol) = _getproperty(edge,Val(s))
_getproperty(edge::EdgeBase,::Val{s}) where s = getfield(edge,s)
_getproperty(edge::EdgeBase,::Val{:nEdges}) = getfield(edge,:n)
_getproperty(edge::EdgeBase,::Val{:verticesOnEdge}) = getfield(edge,:indices).vertices
_getproperty(edge::EdgeBase,::Val{:cellsOnEdge}) = getfield(edge,:indices).cells
_getproperty(edge::EdgeBase,::Val{:xEdge}) = getfield(edge,:position).x
_getproperty(edge::EdgeBase,::Val{:yEdge}) = getfield(edge,:position).y
_getproperty(edge::EdgeBase,::Val{:zEdge}) = getfield(edge,:position).z

struct EdgeVelocityReconstruction{MAX_N_EDGES,TI<:Integer,TF}
    """Number of edges involved in reconstruction of tangential velocity for an edge"""
    nEdges::Vector{TI}
    """IDs of edges involved in reconstruction of tangential velocity for an edge"""
    indices::Vector{VariableLengthStaticVector{MAX_N_EDGES,TI}}
    """Weights used in reconstruction of tangential velocity for an edge"""
    weights::Vector{VariableLengthStaticVector{MAX_N_EDGES,TF}}
end

max_n_edges_vel_reconstruction(::Type{<:EdgeVelocityReconstruction{N}}) where N = N

Base.getproperty(edge::EdgeVelocityReconstruction,s::Symbol) = _getproperty(edge,Val(s))
_getproperty(edge::EdgeVelocityReconstruction,::Val{s}) where s = getfield(edge,s)
_getproperty(edge::EdgeVelocityReconstruction,::Val{:nEdgesOnEdge}) = getfield(edge,:nEdges)
_getproperty(edge::EdgeVelocityReconstruction,::Val{:edgesOnEdge}) = getfield(edge,:indices)
_getproperty(edge::EdgeVelocityReconstruction,::Val{:weightsOnEdge}) = getfield(edge,:weights)

mutable struct EdgeInfo{S,N,TI<:Integer, TF<:Real,Tz<:Number}
    const base::EdgeBase{S,TI,TF,Tz}
    const velRecon::EdgeVelocityReconstruction{N,TI,TF}
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
    normalVectors::TensorsLite.VecMaybe2DxyArray{TF,Tz,1}
    """Weights for cell-centered second derivative, normal to edge, for transport scheme"""
    derivTwo::Array{TF,3}

    function EdgeInfo(edge::EdgeBase{S,TI,TF,Tz},velRecon::EdgeVelocityReconstruction{N}) where {S,N,TI,TF,Tz}
        return new{S,N,TI,TF,Tz}(edge,velRecon)
    end

end

on_a_sphere(::Type{<:EdgeInfo{B}}) where B = B
max_n_edges_vel_reconstruction(::Type{<:EdgeInfo{B,N}}) where {B,N} = N
integer_precision(::Type{<:EdgeInfo{B,N,T}}) where {B,N,T} = T
float_precision(::Type{<:EdgeInfo{B,N,T,TF}}) where {B,N,T,TF} = TF

Base.getproperty(edge::EdgeInfo,s::Symbol) = _getproperty(edge,Val(s))
_getproperty(edge::EdgeInfo,::Val{s}) where s = getfield(edge,s)

for s in (:nEdges,:n,:indices,:position,:verticesOnEdge,:cellsOnEdge,:xEdge,:yEdge,:zEdge)
    @eval _getproperty(edge::EdgeInfo,::Val{$(QuoteNode(s))}) = getproperty(getfield(edge,:base),$(QuoteNode(s)))
end

for s in (:nEdgesOnEdge, :edgesOnEdge, :weightsOnEdge)
    @eval _getproperty(edge::EdgeInfo,::Val{$(QuoteNode(s))}) = getproperty(getfield(edge,:velRecon),$(QuoteNode(s)))
end

for (s,nc) in pairs((longitude=:lonEdge, latitude=:latEdge,
                     indexToID=:indexToEdgeID, dv=:dvEdge, bdyMask=:bdyMaskEdge,
                     dc=:dcEdge, angle=:angleEdge,
                     normalVectors=:edgeNormalVectors, derivTwo=:deriv_two))

    @eval _getproperty(edge::EdgeInfo,::Val{$(QuoteNode(nc))}) = getfield(edge,$(QuoteNode(s)))
end
