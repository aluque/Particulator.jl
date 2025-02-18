#=
Callbacks allow to perform actions during the simulation.
A callback instance must define the methods oncollision (called for each collision)
and onadvance (called for each time the particle is advanced.
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
A combination of several callbacks, called in order.
"""
struct CombinedCallback{T<:Tuple}
    tpl::T
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



    


