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

function init!(mpopl::MultiPopulation, tpl::NamedTuple)
    popl = first(tpl)
    for i in popl.iup[]:popl.n[]
        l = LazyRow(popl.particles, i)
        l.active || continue
        setr!(popl, i)
    end

    init!(mpopl, Base.tail(tpl))
end

function init!(mpopl::MultiPopulation, tpl::Nothing=nothing)
    init!(mpopl, mpopl.index)
end

init!(mpopl::MultiPopulation, tpl::NamedTuple{()}) = nothing


function advance!(mpopl, efield, bfield, tfinal, tracker=VoidCollisionTracker())
    advance_init!(mpopl.index)
    
    # n is the total number particles created during this time-step. We iterate until the number is
    # zero, each time passing only through the particles that have not been updated yet.
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
        
        trem = tfinal - l.t
        while trem > eps(typeof(trem)) && l.active
            tnextcoll = l.s / l.r
            if trem > tnextcoll
                Δt = tnextcoll
                collides = true
            else
                Δt = trem
                collides = false
                l.s -= Δt * l.r
            end
            state = instantiate(l)
            new_state = advance_free(state, efield, bfield, Δt)
            popl.particles[i] = new_state
            
            if collides
                do_one_collision!(mpopl, popl.collisions, new_state, i, tracker)
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
    popl = first(tpl)
    popl.iup[] = 1

    @batch for i in popl.iup[]:popl.n[]
        l = LazyRow(popl.particles, i)
        l.active || continue
        setr!(popl, i)
    end
    
    advance_init!(Base.tail(tpl))
end
advance_init!(tpl::NamedTuple{()}) = nothing
