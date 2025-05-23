# swarm.jl : Tue Feb 18 22:59:01 2025
"""
# Swarm
Investigating swarms under an electric field.

## Running the code

```julia
julia> includet("swarm.jl")     # Needs Revise.jl
julia> Swarm.main()
```
"""
module Swarm


using Particulator

using StaticArrays
using StatsBase
using Accessors

using Random
using Constants: co
using Polyester
using Base.Threads: Atomic
using PyPlot

function main(;n_init_particles=1,
              maxp=1000_000,
              ntarget=10000,
              init_energy=7e6 * co.eV,
              dt=2.5e-11,
              efield = 5e5,
              safety = 1.15,
              z1=0.0,
              z2=200.0,
              seed = 0,
              Kthresh = 1e3 * co.eV,
              zwall=nothing,
              composition = Dict("N2" => co.nair * 0.79,
                                 "O2" => co.nair * 0.21)
              )
    Polyester.reset_threads!()
    Random.seed!(seed)
        
    mc2 = co.electron_mass * co.c^2
    Fdt = co.elementary_charge * efield * dt
    
    # Build collision tables
    ecolls = build_electron_collision_table(composition, Fdt; safety)    
    pcolls = build_positron_collision_table(composition, Kthresh, Fdt; safety)
    γcolls = build_photon_collision_table(composition, Kthresh)
    
    # Initialize particles
    K = init_energy
    v0 = SA[0.0, 10.0, co.c * sqrt(1 - (mc2 / (mc2 + K))^2)]
    
    init_particles = map(1:n_init_particles) do _
        x = SA[0.0, 0.0, 0.0]
        ElectronState(x, v0)
    end

    # Construct populations
    electrons = Population(maxp, init_particles, ecolls, Kthresh)
    photons   = Population(maxp, PhotonState{Float64}[], γcolls, Kthresh)
    positrons = Population(maxp, PositronState{Float64}[], pcolls, Kthresh)

    population_index = Pair{Symbol, Any}[:electron => electrons,
                                         :photon => photons,
                                         :positron => positrons]
    mpopl = MultiPopulation(population_index...)
    init!(mpopl)

    efield = DoubleLayerField(z1, z2, SA[0.0, 0.0, -efield])
    bfield = HomogeneousField(SA[0.0, 0.0, 0.0])
    
    counter = CollisionCounter()
    if !isnothing(zwall)
        callback1 = WallCallback{ElectronState{Float64}}(3, zwall)
        callback2 = WallCallback{PhotonState{Float64}}(3, zwall)
        callback3 = WallCallback{PositronState{Float64}}(3, zwall)

        callback = CombinedCallback((callback1, callback2, callback3))
    else
        callback = VoidCallback();
    end

    tstep = 1e-9
    t = 0.0
    for i in 1:300
        run!(mpopl, efield, bfield, t + tstep, dt, callback; output_dt = nothing)
        n = nactives(electrons)
        if n > ntarget
            roulette!(ntarget / n, electrons)
        end
        t += tstep
        @info n
    end
    
    return NamedTuple(Base.@locals)
end


function build_electron_collision_table(comp, args...; kw...)
    # Tuples with density of scatterers, cross section
    processes = [
        (2 * comp["N2"], RelativisticCoulomb(7)),
        (2 * comp["O2"], RelativisticCoulomb(8)),
        (2 * comp["N2"], SeltzerBerger(7)),
        (2 * comp["O2"], SeltzerBerger(8)),
        [(comp["N2"], orb) for orb in ORBITALS["N2"]]...,
        [(comp["O2"], orb) for orb in ORBITALS["O2"]]...]
                 
    return collision_table_from_processes(processes, Electron, args...; kw...)
end

function build_positron_collision_table(comp, tcut, args...; kw...)
    # Tuples with density of scatterers, cross section
    processes = [(2 * comp["N2"], RelativisticCoulomb(7)),
                 (2 * comp["O2"], RelativisticCoulomb(8)),
                 (2 * comp["N2"], Bhaba(7, tcut)),
                 (2 * comp["O2"], Bhaba(8, tcut)),
                 (2 * comp["N2"], PositronAnihilation(7)),
                 (2 * comp["O2"], PositronAnihilation(8))]
                 
    return collision_table_from_processes(processes, Positron, args...; kw...)
end

function build_photon_collision_table(comp, args...; kw...)
    # Tuples with density of scatterers, cross section
    processes = [(2 * comp["N2"], PhotoElectric(7)),
                 (2 * comp["O2"], PhotoElectric(8)),
                 (2 * comp["N2"], BetheHeitler(7)),
                 (2 * comp["O2"], BetheHeitler(8)),
                 (2 * comp["N2"], Compton(7)),
                 (2 * comp["O2"], Compton(8))]
                     
    return collision_table_from_processes(processes, Photon, 0; kw...)
end


"""
Plot the spectrum of a list of particles `p`.  The spectrum is divided by `n`, which is typically
the number of primaries in an avalanche. `bins` specifies the number of bins in the histogram
(log spaced).
"""
function plotspec(p; n=1, bins=200)
    plt.matplotlib.pyplot.style.use("granada")

    bins = (10 .^ LinRange(log10(1), log10(1e5), bins))
    hst = StatsBase.fit(Histogram, energy.(p) ./ (1e3 * co.eV),
                        weights(getfield.(Particulator.instantiate.(p), :w)),
                        bins)
    hst1 = StatsBase.normalize(hst, mode=:density)
    plt.plot(hst1.edges[1][begin:end-1], hst1.weights ./ n)
    plt.xlabel("energy (keV)")
    plt.ylabel("energy spectrum (particles / keV)")
    return NamedTuple(Base.@locals)
end
end

if abspath(PROGRAM_FILE) == @__FILE__
    Swarm.main()
end

