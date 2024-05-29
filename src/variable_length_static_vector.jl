struct VariableLengthStaticVector{MAX_N,T} <: AbstractVector{T}
    data::NTuple{MAX_N,T}
    length::UInt8
end

function VariableLengthStaticVector(data::NTuple{N,T},length=N) where {N,T}
    return VariableLengthStaticVector{N,T}(data,length)
end

function VariableLengthStaticVector{N_MAX}(data::NTuple{N,T},length=N) where {N_MAX,N,T}
    N > N_MAX && throw(DomainError(data))
    N == N_MAX && return VariableLengthStaticVector(data,length)
    return VariableLengthStaticVector((data...,ntuple(i->zero(T),Val(N_MAX-N))...),length)
end

function VariableLengthStaticVector{N_MAX}(v::VariableLengthStaticVector{N,T}) where{N_MAX,N,T}
    l = length(v)
    N <= N_MAX && return VariableLengthStaticVector{N_MAX}(v.data,l)
    l > N_MAX && throw(DomainError(v,"Input length is grater than the maximum of $N_MAX allowed"))
    return VariableLengthStaticVector{N_MAX}(v.data[Base.OneTo(N_MAX)],l)
end

@inline Base.length(d::VariableLengthStaticVector) = Int(d.length)
@inline Base.size(d::VariableLengthStaticVector) = (length(d),)

@inline function Base.getindex(d::VariableLengthStaticVector,i::Integer)
    @boundscheck checkbounds(d,i)
    data = d.data
    return @inbounds data[i]
end

@inline max_length(::VariableLengthStaticVector{N,T}) where {N,T} = N
@inline max_length(::Type{<:VariableLengthStaticVector{N,T}}) where {N,T} = N

function Base.map(f::F,a::VariableLengthStaticVector{N1,T1}) where {F<:Function,N1,T1}
    Tf = Base.promote_op(f,T1)
    l = length(a)
    d = a.data
    func = i -> (i <= l ? (@inbounds f(d[i])) : zero(Tf))
    return VariableLengthStaticVector{N1,Tf}(ntuple(func,Val{N1}()),a.length)
end

function Base.map(f::F,a::VariableLengthStaticVector{N1,T1},b::VariableLengthStaticVector{N2,T2}) where {F<:Function,N1,T1,N2,T2}
    Tf = Base.promote_op(f,T1,T2)
    Nf = min(N1,N2)
    l = min(length(a),length(b))
    da = a.data
    db = b.data
    func = i -> (i <= l ? (@inbounds f(da[i], db[i])) : (zero(Tf)))
    return VariableLengthStaticVector{Nf,Tf}(ntuple(func,Val{Nf}()),UInt(l))
end
