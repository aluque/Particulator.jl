abstract type CollisionProcess end

struct NullCollision <: CollisionProcess end

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



#
# Collision outcomes
#
# Each collision type must return always the same type of collision to ensure
# type stability.

abstract type AbstractOutcome end

# Nothing happens
struct NullOutcome{PS} <: AbstractOutcome end
    
# The colliding particle dissapears
struct RemoveParticleOutcome{PS} <: AbstractOutcome
    state::PS
end

# The state of the colliding particle changes to p
struct StateChangeOutcome{PS} <: AbstractOutcome
    state::PS
end

# A new particle is created with state state2. The state of
# the colliding particle changes to state1
struct NewParticleOutcome{PS1, PS2} <: AbstractOutcome
    state1::PS1
    state2::PS2
end

# A new particle with state p2 replaces an existing particle with state state1.
# This is e.g. when a photon is absorbed and liberates an electron.
struct ReplaceParticleOutcome{PS1, PS2} <: AbstractOutcome
    state1::PS1
    state2::PS2
end

# Two new particles with states p2, p3 replace an existing particle with state state1.
# This is e.g. a photon creates an e-/e+ pair.
struct ReplaceParticlePairOutcome{PS1, PS2, PS3} <: AbstractOutcome
    state1::PS1
    state2::PS2
    state3::PS3
end

# For Poisson photon generations we may want to generate many particles in a single
# event. To sample them we create instances of this.
struct MultiplePhotonOutcome{PS1, C} <: AbstractOutcome
    # The new state of the electron producing the emission
    nphot::Int
    state1::PS1
    photoemit::C
end


@inline collide(c::NullCollision, k::PS, energy) where PS = NullOutcome{PS}()

"""
Set the r field of particle with index `i` using its energy and the safety factor of the population.
"""
function setr!(popl, i)
    l = LazyRow(popl.particles, i)

    eng = energy(l)
    
    r = ratebound(popl.collisions, eng)
    l.r = r
end

"""
   appply!(population, outcome, i)

Apply a `CollisionOutcome` to a population `population`.  `i` is the index
of the colliding particle, which is needed if it experiences a change.
We delegate to `population` handling the creation of a new particle.
"""
@inline function apply!(mpopl, outcome::NullOutcome{PS}, i) where PS
    popl = get(mpopl, particle_type(PS))
    setr!(popl, i)
    l = LazyRow(popl.particles, i)
    l.s = nextcoll()
end

@inline function apply!(mpopl, outcome::StateChangeOutcome{PS}, i) where PS
    popl = get(mpopl, particle_type(PS))
    popl.particles[i] = outcome.state
    setr!(popl, i)
end

@inline function apply!(mpopl, outcome::NewParticleOutcome{PS1, PS2}, i) where {PS1, PS2}
    popl1 = get(mpopl, particle_type(PS1))
    popl1.particles[i] = outcome.state1
    setr!(popl1, i)

    popl2 = get(mpopl, particle_type(PS2))
    j = add_particle!(popl2, outcome.state2)

    j > 0 && setr!(popl2, j)
end

@inline function apply!(mpopl, outcome::RemoveParticleOutcome{PS}, i) where PS
    popl = get(mpopl, particle_type(PS))
    remove_particle!(popl, i)
end

@inline function apply!(mpopl, outcome::ReplaceParticleOutcome{PS1, PS2}, i) where {PS1, PS2}
    popl1 = get(mpopl, particle_type(PS1))
    remove_particle!(popl1, i)

    popl2 = get(mpopl, particle_type(PS2))
    j = add_particle!(popl2, outcome.state2)
    j > 0 && setr!(popl2, j)

end

@inline function apply!(mpopl, outcome::ReplaceParticlePairOutcome{PS1, PS2, PS3}, i) where {PS1, PS2, PS3}
    popl1 = get(mpopl, particle_type(PS1))
    remove_particle!(popl1, i)

    popl2 = get(mpopl, particle_type(PS2))
    i = add_particle!(popl2, outcome.state2)    
    i > 0 && setr!(popl2, i)
    
    popl3 = get(mpopl, particle_type(PS3))
    j = add_particle!(popl3, outcome.state3)    
    j > 0 && setr!(popl3, j)
end


#
# With collision tracker we can define a callback executed whenever a collision
# happens.
#
abstract type AbstractCollisionTracker end
struct VoidCollisionTracker <: AbstractCollisionTracker end

track(::AbstractCollisionTracker, collision, outcome) = nothing


indweight(colls::CollisionTable, E) = indweight(colls.energy, E)

"""
Sample one (possibly null) collision.
"""    
@generated function do_one_collision!(mpopl, colls::AbstractCollisionTable{T, C}, state, i, tracker) where {T, C}
    L = fieldcount(C)

    out = quote
        $(Expr(:meta, :inline))

        eng = energy(state)
        pre = presample(colls, state, eng)
        ξ = rand(T) * state.r

        # Check that we are not underestimating rate bound.
        # rr = zero(state.r)
        # for i in 1:$L
        #     rr += rate(colls, i, pre)
        # end
        # @assert rr <= 1.00001 * state.r "Rate upper bound underestimated ($rr vs $(state.r) energy=$(eng/co.eV) eV $(state)): increase safety factor"
        
        #k, w = indweight(colls, energy)
    end
    
    for j in 1:L
        push!(out.args,
              quote
                  ν = rate(colls, $j, pre)
                  if ν > ξ                      
                      outcome = collide(colls.proc[$j], state, eng)
                      track(tracker, colls.proc[$j], outcome)
                      apply!(mpopl, outcome, i)
                      return
                  else
                      ξ -= ν
                  end
              end)
    end
    
    push!(out.args,
          quote
              # Even after exhausting all processes ξ should remain positive if we
              # correctly bounded the max. rate.
              @assert ξ >= 0 "Rate upper bound underestimated: increase safety factor"

              
              # This is a null collision. Only takes effect if we somehow want to track them or
              # do something weird. Otherwise all this should fold to a nop.
              outcome = collide(NullCollision(), state, eng)
              track(tracker, NullCollision(), outcome)
              apply!(mpopl, outcome, i)

              return nothing
          end)
          
    return out
end
