abstract type AbstractCollisionTable{T, C <: Tuple}; end

"""
Struct with a set of collisions evaluated at the same energy grid.

The rate bound can either be a number or another struct. It is used in the form

ratebound(::CollisionTable, energy).

If rate bound is a number the  assumption is that the total collision rate is independent of
energy (because the space is filled with NullCollisions).

The rate bound at a given energy must be an upper limit of the energy-dependent rate during a timestep.
"""
@kwdef struct CollisionTable{T, C <: Tuple, V <: AbstractRange{T},
                             A <: AbstractMatrix{T}, TR <: Union{Vector{T}, T}} <: AbstractCollisionTable{T, C}
    "Each type of collision"
    proc::C

    "Energy grid"
    energy::V

    "Array with collision rates for each process."
    rate::A

    "Upper limit of the total collision rate during a time-step"
    ratebound::TR
end

Base.length(c::CollisionTable) = length(c.proc)

ratebound(c::CollisionTable, eng, pre=nothing) = ratebound(c.ratebound, c, eng, pre)
ratebound(x::Number, c, eng, pre=nothing) = x

function ratebound(v::Vector, c, eng, pre::Nothing)
    pre = presample(c, nothing, eng)
    ratebound(v, c, eng, pre)
end

function ratebound(v::Vector, c, eng, pre)
    k, w = pre
    return w * v[k] + (1 - w) * v[k + 1]
end

maxenergy(c::CollisionTable) = last(c.energy)

# Checking for collisions involves many tests but some computations are
# common to all of them. Here we store them in a generic way. energy is passed
# as an optimization.
presample(c::CollisionTable, state, energy) = indweight(c.energy, energy)

# Obtain probability rate of process j
function rate(c::CollisionTable, j, pre)
    k, w = pre

    return w * c.rate[j, k] + (1 - w) * c.rate[j, k + 1]
end


"""
A collision table based on Chebyshev expansions in a set of binary intervals.
"""
@kwdef struct ChebyshevCollisionTable{T, N, C <: Tuple, A <: AbstractArray, B} <: AbstractCollisionTable{T, C}
    "Each type of collision"
    proc::C

    "Binary divisions"
    b::BinaryIntervals{T}

    "Chebyshev expansions for each interval for each process"
    rate::A

    "Chebyshev expansions for the upper limit of the total collision rate during a time-step"
    ratebound::B
end


Base.length(c::ChebyshevCollisionTable) = length(c.proc)

maxenergy(c::ChebyshevCollisionTable) = c.b.xmax

presample(c::ChebyshevCollisionTable{T, N}, state, energy) where {T, N} = precheb(energy, c.b, Val{N}())

# Obtain probability rate of process j
function rate(c::ChebyshevCollisionTable{T, N}, j, pre) where {T, N}
    (i, t) = pre

    # i = interval
    # j = process
    # k = order
    @inbounds return sum(ntuple(k -> @inbounds(c.rate[k, j, i + 1] * t[k]), Val(N)))    
end

function ratebound(c::ChebyshevCollisionTable, eng, pre::Nothing=nothing)
    pre = presample(c, nothing, eng)
    ratebound(c, eng, pre)
end

function ratebound(c::ChebyshevCollisionTable{T, N}, eng, pre) where {T, N}
    (i, t) = pre

    # i = interval
    # j = process
    # k = order
    return sum(ntuple(k -> c.ratebound[k, i + 1] * t[k], Val(N)))    
end
