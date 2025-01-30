""" 
Sample points from the unitary sphere. 
"""
function randsphere()
    ϕ = 2π * rand()
    sinϕ, cosϕ = sincos(ϕ)
    
    u = 2 * rand() - 1
    v = sqrt(1 - u^2)

    @SVector [v * cosϕ, v * sinϕ, u]
end

"""
Sample the next collision time for a particle with unit collision rate.
"""
nextcoll() = -log(rand())


"""
Return the index and weight of the item before `x` in a range `r`.
"""
function indweight(r::LinRange, x)
    #i = searchsortedlast(r, x)
    #@show r i
    i = Int(fld(x - first(r), step(r))) + 1
    #@assert ip == i "$ip != $i"
    
    w = (r[i + 1] - x) / step(r)
    @assert 0 <= w <= 1
    return (i, w)
end



""" 
Find a vector with the same norm as `u` but deviated with inclination `θ` and azimuth
`ϕ` with respect to u.  For performance, the passed parameter is `cosθ` instead of `θ`.
"""
@inline @fastmath function turn(u, cosθ, ϕ)
    μ = u / norm(u)
    sinϕ, cosϕ = sincos(ϕ)
    
    sinθ = sqrt(1 - cosθ^2)

    s = sqrt(1 - μ[3]^2)

    b(μ1, μ2, μ3) = μ1 * μ3 * cosϕ + μ2 * sinϕ 

    μ1x = sinθ * b(μ[1], -μ[2], μ[3]) / s + μ[1] * cosθ
    μ1y = sinθ * b(μ[2],  μ[1], μ[3]) / s + μ[2] * cosθ
    μ1z = -s * sinθ * cosϕ + μ[3] * cosθ

    v = norm(u) * SA[μ1x, μ1y, μ1z]

    return v
end


struct LogLinRange{T} <: AbstractRange{T}
    L::LinRange{T, Int}
    x0::T
end

"""
Construct a lin-log scale with the form x = exp(L) - x0, where L is a linear range. The range is
constructed such that for small x the intervals in x are about `Δx`. The max x value is about `xmax` and
there are `N` points.
"""
function LogLinRange(L1, L2, N) 
    L = LinRange(L1, L2, N)
    x0 = exp(first(L))

    return LogLinRange(L, x0)
end

Base.getindex(r::LogLinRange, i::Integer) = exp(r.L[i]) - r.x0
Base.length(r::LogLinRange) = length(r.L)
Base.lastindex(r::LogLinRange) = lastindex(r.L)
Base.firstindex(r::LogLinRange) = firstindex(r.L)
Base.isempty(r::LogLinRange) = length(r) == 0
Base.step(::LogLinRange) = error("LogLinRange does not have a defined step")
Base.first(r::LogLinRange) = exp(first(r.L)) - r.x0
Base.last(r::LogLinRange) = exp(last(r.L)) - r.x0
Base.argmin(r::LogLinRange) = firstindex(r)
Base.argmax(r::LogLinRange) = lastindex(r)
Base.axes(r::LogLinRange) = axes(r.L)

# Prevent specializations for LinRanges
for op in [:+, :-, :*, :/]
    @eval begin
        Base.broadcasted(s::Base.Broadcast.DefaultArrayStyle{1}, ::typeof($op), x::Number, r::LogLinRange) = Base.Broadcast.Broadcasted(s, $op, (x, r))
        Base.broadcasted(s::Base.Broadcast.DefaultArrayStyle{1}, ::typeof($op), r::LogLinRange, x::Number) = Base.Broadcast.Broadcasted(s, $op, (r, x))
    end
end

function Base.iterate(r::LogLinRange, i::Integer=zero(length(r)))
    @inline
    i += oneunit(i)
    length(r) < i && return nothing
    getindex(r, i), i
end

for s in [:searchsortedfirst, :searchsortedlast, :searchsorted]
    @eval begin
        Base.$s(r::LogLinRange, x;
           lt=isless, by=identity, rev::Union{Bool,Nothing}=nothing, order::Base.Ordering=Base.Forward) =
               Base.$s(r.L, log(x + r.x0); lt, by, rev, order)
    end
end

            
function Base.show(io::IO, r::LogLinRange)
    print(io, "LogLinRange(", repr(r.L), ", ", repr(r.x0), ")")
end


function indweight(r::LogLinRange, x)
    l = log(x + r.x0)
    
    i = Int(fld(l - first(r.L), step(r.L))) + 1
    #@assert ip == i "$ip != $i"

    w = (exp(r.L[i + 1]) - r.x0 - x) / (exp(r.L[i + 1]) - exp(r.L[i]))
    #@assert 0 <= w <= 1
    return (i, w)
end


#
# Reading of Geant4 2d vectors
#
struct RawG4Physics2DVector{T}
    k::Int
    nx::Int
    ny::Int

    x::Vector{T}
    y::Vector{T}
    value::Matrix{T}
end

function RawG4Physics2DVector{T}(io::IO) where T
    (k, nx, ny) = parse.(Int, split(readline(io)))
    x = parse.(T, split(readline(io)))
    y = parse.(T, split(readline(io)))

    @assert length(x) == nx
    @assert length(y) == ny
    
    value = readdlm(io, T)

    @assert size(value) == (ny, nx)

    return RawG4Physics2DVector{T}(k, nx, ny, x, y, value)
end
