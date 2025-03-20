#=
Callbacks allow to perform actions during the simulation.
A callback instance must define the methods oncollision (called for each collision),
onadvance (called for each time the particle is advanced, and onstep (called after
each completed timestep).
=#

abstract type AbstractCallback end
struct VoidCallback <: AbstractCallback end

"""
Callback method called whenever there is a a `collision` happens with `outcome`.
If another `outcome` is returned, it will be used instead of the original.
"""
oncollision(::AbstractCallback, collision, outcome, t=nothing) = outcome

"""
Callback method called whenever a particle is advanced. `old_state`is the particle state before
the advance and `new_state` is the state after the advance. The callback must return
a particle state that will be used instead of `new_state` if it differs from it.
"""
onadvance(::AbstractCallback, old_state, new_state, t=nothing) = new_state

"""
Callback method called after each timestep. `mpopl` is a `MultiPopulation` containing
the current state of all particles from each population.
"""
onstep(::AbstractCallback, mpopl, t=nothing) = nothing

"""
Callback method called at each output timestep. `mpopl` is a `MultiPopulation` containing
the current state of all particles from each population.
"""
onoutput(::AbstractCallback, mpopl, t=nothing, i=nothing) = nothing


"""
A combination of several callbacks, called in order.
"""
struct CombinedCallback{T<:Tuple} <: AbstractCallback
    tpl::T
end

function CombinedCallback(c::CombinedCallback, other::Tuple)
    CombinedCallback((c.tpl..., other...))
end

function CombinedCallback(c::VoidCallback, other::Tuple)
    CombinedCallback(other)
end

function CombinedCallback(c::AbstractCallback, other::Tuple)
    CombinedCallback((c, other...))
end

Base.getindex(c::CombinedCallback, i::Int) = c.tpl[i::Int]

function oncollision(c::CombinedCallback, args...)
    _oncollision(c.tpl, args...)
end

function _oncollision(tpl::Tuple, collision, outcome, t=nothing)
    f = first(tpl)
    outcome = oncollision(f, collision, outcome, t)
    return _oncollision(Base.tail(tpl), collision, outcome, t)
end

_oncollision(tpl::Tuple{}, collision, outcome, t=nothing) = outcome

function onadvance(c::CombinedCallback, args...)
    _onadvance(c.tpl, args...)
end

function _onadvance(tpl::Tuple, old_state, new_state, t=nothing)
    f = first(tpl)
    new_state = onadvance(f, old_state, new_state, t)
    return _onadvance(Base.tail(tpl), old_state, new_state, t)
end

_onadvance(tpl::Tuple{}, old_state, new_state, t=nothing) = new_state

function onoutput(c::CombinedCallback, args...)
    _onoutput(c.tpl, args...)
end

function _onoutput(tpl::Tuple, mpopl, t=nothing, i=nothing)
    f = first(tpl)
    onoutput(f, mpopl, t, i)
    _onoutput(Base.tail(tpl), mpopl, t, i)
    return nothing
end

_onoutput(tpl::Tuple{}, mpopl, t=nothing, i=nothing) = nothing



#
# Some useful callbacks
#
"""
A callback that counts collisions and organizes them according to collision type.
"""
struct CollisionCounter <: AbstractCallback
    d::Dict{Type, Atomic{Int}}
    CollisionCounter() = new(Dict{Type, Int}())
end

function oncollision(cc::CollisionCounter, collision, outcome, t)
    @info "Collision at t=$t"
    s = typeof(collision)
    if haskey(cc.d, s)
        cc.d[s][] += 1
    else
        cc.d[s] = Atomic{Int}(1)
    end
    return outcome
end

function Base.show(io::IO, cc::CollisionCounter)
    println(io, "CollisionCounter with values:")
    tot = 0
    for (key, value) in cc.d
        println(io, styled"{magenta:$key}:   {red:$(value[])}")
        tot += value[]
    end
    println(io, styled"Total:   {green:$(tot)}")
end

"""
A callback to trace particles as they cross a wall.
"""
struct WallCallback{P, T, L} <: AbstractCallback
    coord::Int
    v::T
    drop::Bool

    accum::L
    lck::ReentrantLock
    
    """
    Initialize a WallCallback to count particles with states of type `P`. The wall is perpendicular
    to the axis `coord` and located at `v` (e.g. z=1 if `coord=3`, `v=1`). `drop` indicates whether
    the particle must be dropped after recording it (defaults to `true`).
    """
    function WallCallback{P}(coord, v, drop=true) where P
        accum = P[]
        sizehint!(accum, 10000)
        new{P, typeof(v), typeof(accum)}(coord, v, drop, accum, ReentrantLock())
    end
end


function Particulator.onadvance(wcb::WallCallback{P}, old_state::P, new_state::P, t) where P
    (;coord, v, accum, drop, lck) = wcb
    if old_state.x[coord] < v < new_state.x[coord]
        # interpolate
        w = (v - old_state.x[coord]) / (new_state.x[coord] - old_state.x[coord])
        #
        let mid_state = lincomb(new_state, old_state, w)
            lock(lck) do
                push!(accum, mid_state)
            end
        end
        if drop
            @reset new_state.active = false
        end
    end

    return new_state
end


"""
A callback to save particle numbers.
"""
struct ParticleCountCallback <: AbstractCallback
    p::Vector{Symbol}
    counts::Vector{Vector{Float64}}
    
    function ParticleCountCallback(particles)
        counts = Vector{Vector{Float64}}()
        p = copy(particles)
        return new(p, counts)
    end
end


function onoutput(pc::ParticleCountCallback, mpopl, t, i)
    push!(pc.counts, [t; (nparticles(get(mpopl, ParticleType{s})) for s in pc.p)...])           
end
