# Managing of a populations of homogeneous particles.

# Population of a given particle type
# A note about iup: this is introduced in order to advance particles created during a timestep.
# When we start the advancing, it has to be 1; then after going to all particles it's set to the
# index of the last particle that we have already advanced.
struct Population{PS <: ParticleState, T, C, K, I}
    "Number of particles"
    n::Atomic{Int}

    iup::Atomic{Int}
    
    "Vector with the particles"
    particles::StructArray{PS, 1, K, I}

    "Collision table"
    collisions::C

    "Energy cut of this population"
    energy_cut::T
end


# Ensure that something is not a LazyRow
instantiate(ps::ParticleState) = ps
instantiate(ps::LazyRow{P}) where P <: ParticleState = getfield(ps, 1)[getfield(ps, 2)]


"""
Initialize a population with space to contain `max_particles`,
starting with a vector of super-particle states.
The collision table must be also given in `collisions`.
"""
function Population(max_particles::Int, inparticles::Vector{PS},
                    collisions, energy_cut=0.0) where PS <: ParticleState
    init_particles = length(inparticles)
    particles = StructVector{PS}(undef, max_particles)

    for i in 1:init_particles
        particles[i] = inparticles[i]
    end

    return Population(Atomic{Int}(init_particles), Atomic{Int}(1), particles, collisions, energy_cut)
end

struct ParticleIterator{P <: Population}
    p::P
end
function Base.iterate(iter::ParticleIterator, state=1)
    iter.p.n[] >= state ? (LazyRow(iter.p.particles, state), state + 1) : nothing
end

struct ActiveParticleIterator{P <: Population}
    p::P
end

function Base.iterate(iter::ActiveParticleIterator, i=1)
    while i <= iter.p.n[]
        if iter.p.particles.active[i]
            return (LazyRow(iter.p.particles, i), i + 1)
        end
        i += 1
    end
    return nothing
end

Base.length(iter::ActiveParticleIterator) = nactives(iter.p)

function Base.empty!(popl::Population)
    popl.n[] = 0
end
    
isactive(popl::Population, i::Int) = popl.particles.active[i]

particle_type(popl::Population{PS}) where {PS} = particle_type(PS)

"Return number of particles (active or not)."
nparticles(popl::Population) = popl.n[]

"Iterate over all particles."
eachparticle(popl::Population) = LazyRows(view(popl.particles, 1:popl.n[]))

"Iterate over active particles."
eachactive(popl::Population) = ActiveParticleIterator(popl)

"""
Return number of active particles in the population.
"""
function nactives(popl::Population)
    n = 0
    for i in 1:popl.n[]
        if LazyRow(popl.particles, i).active
            n += 1
        end
    end
    return n
end


"""
Add a particle to the population `popl` with super-state super_state.
"""
function add_particle!(popl::Population, state::ParticleState)
    (;n, particles, collisions) = popl
    (kinenergy(state) <= popl.energy_cut) && return -1
    
    @assert n[] < length(particles) "Maximum number of particles reached"

    nprev = atomic_add!(n, 1)
    
    particles[nprev + 1] = state
    return nprev + 1
end


""" 
Remove particle `i` from the population `popl` by setting its `active` flag
to false
"""
function remove_particle!(popl::Population, i::Integer)
    LazyRow(popl.particles, i).active = false
end


_particle_state_type_param(::Type{<:ParticleState{T}}) where T = T

"""
Compute the total weight of a population `popl`.
"""
function weight(popl::Population{PS}) where PS
    T = _particle_state_type_param(PS)
    weight = zero(T)
    for p in eachactive(popl)
        if p.active
            weight += p.w
        end
    end

    return weight
end


"""
Compute the number of active particles in a population `pop`.
"""
actives(popl) = count(p -> p.active, eachparticle(popl))


"""
Compute the mean energy of a population `popl`.
"""
function meanenergy(popl::Population{PS}) where PS
    T = _particle_state_type_param(PS)

    tot = zero(T)
    totw = zero(T)
    
    nparts = 0
    for p in eachparticle(popl)
        if p.active
            totw += p.w
            tot += p.w * kinenergy(instantiate(p))
        end
    end
    tot / totw
end


"""
Compute the highest energy of a population `popl`.
"""
function maxenergy(popl::Population{PS}) where PS
    maximum(p -> kinenergy(instantiate(p)), eachparticle(popl))
end


"""
Compute the centroid and deviation of a population `popl`.
"""
function spread(popl::Population{PS}) where PS
    T = _particle_state_type_param(PS)

    xmean = @SVector zeros(T, 3)
    x2mean = zero(T)
    tot = zero(T)
    
    nparts = 0
    for p in eachparticle(popl)
        if p.active
            xmean += p.w * p.x
            x2mean += p.w * dot(p.x, p.x)
            tot += p.w
        end
    end
    xmean /= tot
    x2mean /= tot
    return (xmean, sqrt(abs(x2mean - dot(xmean, xmean))))
end


"""
Compute the variances in location of a population `popl`.
"""
function posvar(popl::Population{PS}) where PS
    T = _particle_state_type_param(PS)

    xsum = @SVector zeros(T, 3)
    x2sum = @SVector zeros(T, 3)
    tot = zero(T)
    
    nparts = 0
    for p in eachparticle(popl)
        if p.active
            xsum += p.w * p.x
            x2sum += p.w .* p.x .* p.x
            tot += p.w
        end
    end
    x2sum ./= tot
    xsum ./= tot

    return x2sum .- xsum .^ 2
end

"""
Reorders the particles in the population `popl` to have all active particle 
at the initial positions in the list.
"""
function repack!(popl::Population)
    # last active
    l = nparticles(popl)

    l == 0 && return

    while l > 0 && !isactive(popl, l)
        l -= 1
    end

    if l == 0
        popl.n[] = 0
        return
    end
    
    i = 1
    while i <= l
        if !isactive(popl, i)
            popl.particles[i] = popl.particles[l]
            l -= 1
            while !isactive(popl, l)
                l -= 1
            end
        end
        i += 1
    end

    popl.n[] = l

    return
end


""" 
Shuffle a population using a permutation sampled from a uniform distribution
of permutations.
"""
function shuffle!(popl::Population)
    for i in nparticles(popl):-1:2
        j = rand(1:i)
        popl.particles[i], popl.particles[j] = popl.particles[j], popl.particles[i]
    end        
end

function droplow!(popl, thres=0.0)
    thres = thres == 0 ? popl.energy_cut : thres
    
    prt = popl.particles
    for i in 1:popl.n[]
        if prt[i].active && kinenergy(prt[i]) < thres
            remove_particle!(popl, i)
        end
    end

    repack!(popl)
end


"""
Apply Russian roulette to the population. `f` is a function that receives the
particle energy and returns a probabilty that the particle be retained.
"""
function roulette!(f::Function, popl::Population)
    n = 0
    for i in 1:popl.n[]        
        l = LazyRow(popl.particles, i)
        l.active || continue
        eng = kinenergy(instantiate(l))
        p = f(eng)
        if rand() < p
            l.w /= p
        else
            remove_particle!(popl, i)
        end
    end
end

"""
Apply Russian roulette with a constant probability of dropping a particle.
"""
roulette!(f::Number, popl::Population) = roulette!(_ -> f, popl)


"""
Apply splitting to the population. `f` is a function that receives the
particle energy and returns the mean number of new particles that will be split from it.
"""
function split!(f::Function, popl::Population)
    n = 0
    for i in 1:popl.n[]        
        l = LazyRow(popl.particles, i)
        l.active || continue
        p = f(instantiate(l))
        l.w /= (1 + p)

        state = instantiate(l)
        eng = kinenergy(state)
        for i in 1:rand(Poisson(p))
            add_particle!(popl, state)
        end
    end
end

"""
Apply Russian roulette with a constant probability of dropping a particle.
"""
split!(f::Number, popl::Population) = split!(_ -> f, popl)

