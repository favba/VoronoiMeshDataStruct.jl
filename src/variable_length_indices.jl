
struct VariableLengthIndices{MAX_N,T<:Integer} <: AbstractVector{T}
    data::NTuple{MAX_N,T}

    function VariableLengthIndices(data::NTuple{N,T}) where {T,N}
        i = 1
        val = data[i]
        while val > zero(T) && i < N
            i += 1
            val = data[i]
        end
        if i != N
            all(x->(x==zero(T)), data[i:N]) || throw(DomainError(data,"Invalid indices tuple"))
        end
        return new{N,T}(data)
    end

end

function VariableLengthIndices{N_MAX}(data::NTuple{N,T}) where {N_MAX,N,T}
    N > N_MAX && throw(DomainError(data))
    N == N_MAX && return VariableLengthIndices(data)
    return VariableLengthIndices((data...,ntuple(i->zero(T),Val(N_MAX-N))...))
end

function VariableLengthIndices{N_MAX}(v::VariableLengthIndices{N,T}) where{N_MAX,N,T}
    N <= N_MAX && return VariableLengthIndices{N_MAX}(v.data)
    l = length(v)
    l > N_MAX && throw(DomainError(v,"Input length is grater than the maximum of $N_MAX allowed"))
    return VariableLengthIndices{N_MAX}(v.data[Base.OneTo(N_MAX)])
end

@inline function Base.size(d::VariableLengthIndices{N,T}) where {T,N}
    data = d.data
    l = N
    @inbounds v = data[l]
    while v == zero(T) && l > 0
        l-=1
        @inbounds v = data[l]
    end
    return (l,)
end

@inline function Base.getindex(d::VariableLengthIndices,i::Integer)
    @boundscheck checkbounds(d,i)
    data = d.data
    return @inbounds data[i]
end

@inline Base.iterate(d::VariableLengthIndices) = (d.data[1],2)

@inline function Base.iterate(d::VariableLengthIndices{N,T},state) where {T,N}
    if state > N
        return nothing
    else
        @inbounds val = d.data[state]
        return val == zero(T) ? nothing : (val, state+1)
    end
end

@inline max_length(::VariableLengthIndices{N,T}) where {N,T} = N
@inline max_length(::Type{VariableLengthIndices{N,T}}) where {N,T} = N
