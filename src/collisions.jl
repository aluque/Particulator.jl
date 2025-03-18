abstract type CollisionProcess end

struct NullCollision <: CollisionProcess end

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
    if eng < popl.energy_cut
        l.r = 0
    else
        r = ratebound(popl.collisions, eng)
        l.r = r
    end
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



indweight(colls::CollisionTable, E) = indweight(colls.energy, E)

"""
Sample one (possibly null) collision.
"""    
@generated function do_one_collision!(mpopl, popl, colls::AbstractCollisionTable{T, C},
                                      state, i, t, callback) where {T, C}
    L = fieldcount(C)
    
    out = quote
        $(Expr(:meta, :inline))
        state.r == 0 && return nothing
        
        eng = energy(state)
        eng >= popl.energy_cut || return nothing
        
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
                      outcome = oncollision(callback, colls.proc[$j], outcome, t)
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
              oncollision(callback, NullCollision(), outcome)
              apply!(mpopl, outcome, i)

              return nothing
          end)
          
    return out
end
