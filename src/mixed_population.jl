"""
A mixed population, composed of multiple particle types.
"""
struct MultiPopulation{TP <: NamedTuple}
    # index is a named tuple that connects particle types with their population
    index::TP

    function MultiPopulation(d...)
        index = NamedTuple(d)
        new{typeof(index)}(index)
    end
end


Base.get(mp::MultiPopulation, pt::Type{ParticleType{S}}) where S = getfield(mp.index, S)
Base.map(f, mp::MultiPopulation) = map(f, mp.index)
Base.pairs(mp::MultiPopulation) = pairs(mp.index)

function advance!(mpopl, efield, bfield, Δt, tracker=VoidCollisionTracker())
    advance1!(mpopl.index, mpopl, efield, bfield, Δt, tracker)
end

"""
Advance the particles in the population performing, if needed, intermediate
collisions.
"""
function advance1!(tpl::NamedTuple, mpopl, efield, bfield, Δt, tracker=VoidCollisionTracker())
    popl = first(tpl)
    @batch for i in 1:popl.n[]
        l = LazyRow(popl.particles, i)
        l.active || continue
        
        trem = Δt
        while trem > 0
            tnextcoll = l.s / maxrate(popl.collisions)
            if trem > tnextcoll
                t = tnextcoll
                collides = true
            else
                t = trem
                collides = false
                l.s -= t * maxrate(popl.collisions)
            end
            state = instantiate(l)
            new_state = advance_free(state, efield, bfield, t)
            popl.particles[i] = new_state
            
            if collides
                do_one_collision!(mpopl, popl.collisions, new_state, i, tracker)
                l.s = nextcoll()
            end
            trem -= t
        end
    end

    advance1!(Base.tail(tpl), mpopl, efield, bfield, Δt, tracker)
end

advance1!(tpl::NamedTuple{()}, mpopl, efield, bfield, Δt, tracker=VoidCollisionTracker()) = nothing


"""
Check for the particles that are set to collide and then perform a
(possibly null) collision
"""
function collisions!(mpopl, Δt, tracker=VoidCollisionTracker())
    # WARN: Possibly type-unstable
    for (sym, popl) in pairs(mpopl)
        @batch for i in 1:popl.n[]
            l = LazyRow(popl.particles, i)
            l.active || continue

            l.s -= Δt * maxrate(popl.collisions)
            if l.s <= 0
                state = instantiate(l)
                do_one_collision!(mpopl, popl.collisions, state, i, tracker)
                l.s = nextcoll()
            end
        end
        #@info "$c/$m = $(c/m) collision fraction"
    end
end
