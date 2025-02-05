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

function advance!(mpopl, efield, bfield, tfinal, tracker=VoidCollisionTracker())
    advance_init!(mpopl.index)
    
    local n = 1
    while n > 0
        n = advance1!(mpopl.index, mpopl, efield, bfield, tfinal, tracker)
    end
end

"""
Advance the particles in the population performing, if needed, intermediate
collisions.
"""
function advance1!(tpl::NamedTuple, mpopl, efield, bfield, tfinal, tracker=VoidCollisionTracker())
    popl = first(tpl)
    (;collisions, iup) = popl
    ilast = popl.n[]

    @batch for i in iup[]:ilast
        l = LazyRow(popl.particles, i)
        l.active || continue
        
        eng = energy(instantiate(l))
        totrate = totalrate(collisions, eng)
        
        trem = tfinal - l.t
        while trem > eps(typeof(trem))
            tnextcoll = l.s / totrate
            if trem > tnextcoll
                Δt = tnextcoll
                collides = true
            else
                Δt = trem
                collides = false
                l.s -= Δt * totrate
            end
            state = instantiate(l)
            new_state = advance_free(state, efield, bfield, Δt)
            popl.particles[i] = new_state
            
            if collides
                do_one_collision!(mpopl, popl.collisions, new_state, i, tracker)
                l.s = nextcoll()
            end
            trem -= Δt
        end
    end
    n = ilast - iup[] + 1
    iup[] = ilast + 1
    
    return n + advance1!(Base.tail(tpl), mpopl, efield, bfield, tfinal, tracker)
end

advance1!(tpl::NamedTuple{()}, mpopl, efield, bfield, tfinal, tracker=VoidCollisionTracker()) = 0

function advance_init!(tpl::NamedTuple)
    first(tpl).iup[] = 1
    advance_init!(Base.tail(tpl))
end
advance_init!(tpl::NamedTuple{()}) = nothing


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

            eng = energy(instantiate(l))
            totrate = totalrate(colls, eng)

            l.s -= Δt * totrate
            if l.s <= 0
                state = instantiate(l)
                do_one_collision!(mpopl, popl.collisions, state, i, tracker)
                l.s = nextcoll()
            end
        end
        #@info "$c/$m = $(c/m) collision fraction"
    end
end
