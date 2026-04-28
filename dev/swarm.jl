module Swarm

using StaticArrays
using StatsBase

using Random
using Constants: co
using Polyester
using Base.Threads: Atomic
using DataFrames

using Particulator

using Particulator
import Particulator: ParticleType

function main(;
              n_init_particles=1,
              ntarget=10000,
              maxp=max(1000_000, 10 * n_init_particles),
              init_energy=1e6 * co.eV,
              dt=2.5e-11,
              efield = 3.5e5,
              safety = 1.15,
              tfinal = 1e-8,
              z1=0.0,
              z2=2000.0,
              seed = 0,
              Kthresh = 1e3 * co.eV,
              output_dt = nothing,
              zwall = nothing,
              density_factor=1.0,
              composition = Dict("N2" => co.nair * 0.79 * density_factor,
                                 "O2" => co.nair * 0.21 * density_factor))
    Polyester.reset_threads!()
    Random.seed!(seed)
        
    mc2 = co.electron_mass * co.c^2
    Fdt = co.elementary_charge * efield * dt
    
    # Build collision tables
    ecolls = build_electron_collision_table(composition, Kthresh, Fdt; safety)    
    pcolls = build_positron_collision_table(composition, Kthresh, Fdt; safety)
    γcolls = build_photon_collision_table(composition)
    
    # Initialize particles
    K = init_energy
    pnorm = Particulator.momentum_norm_from_kin(Electron, K)
    p0 = SA[0.0, 1e-6 * pnorm, pnorm]
    
    init_particles = map(1:n_init_particles) do _
        x = SA[0.0, 0.0, 0.0]
        ElectronState(x, p0)
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

    #efield = ConfinedDoubleLayerField(1e3, 1e3, 500.0, -efield)    
    efield = DoubleLayerField(z1, z2, SA[0.0, 0.0, -efield])
    bfield = HomogeneousField(SA[0.0, 0.0, 0.0])
    emfield = ElectromagneticField(efield, bfield)


    # I from M.J. Berger et al. Icru report 37. Journal of the International
    # Commission on Radiation Units and Measurements, os19(2):, dec 1984
    I = ((composition["N2"] * 82 + composition["O2"] * 95) * co.eV /
        (composition["N2"] + composition["O2"]))

    # Electron density
    nel = composition["N2"] * 14 + composition["O2"] * 16
    loss1 = ContinuumLoss(nel, I, Kthresh)
    loss = ChebContinuumLoss(loss1, 1e9 * co.eV, 3)

    loss = RestrictedForcing{PositronState}(loss)
    
    pusher = RK2Pusher(CombinedForcing(emfield, loss))
    #pusher = RK2Pusher(emfield)
    
    
    counter = CollisionCounter()
    if !isnothing(zwall)
        callback1 = WallCallback{ElectronState{Float64}}(3, zwall)
        callback2 = WallCallback{PhotonState{Float64}}(3, zwall)
        callback3 = WallCallback{PositronState{Float64}}(3, zwall)

        callback = CombinedCallback((callback1, callback2, callback3))
    else
        callback = VoidCallback();
    end

    particle_count = ParticleCountCallback([:electron, :photon, :positron])
    thinning = RouletteCallback(ntarget)

    callback = CombinedCallback(callback, (particle_count, thinning))
    #callback = CombinedCallback(callback, (particle_count,))
    
    run!(mpopl, pusher, tfinal, dt, callback; output_dt)
    
    counts = DataFrame([s => [v[i] for v in particle_count.counts] for (i, s) in
                            enumerate([:t, :electron, :photon, :positron])])
    
    return NamedTuple(Base.@locals)
end


function build_electron_collision_table(comp, tcut, args...; kw...)
    # Tuples with density of scatterers, cross section
    processes = [
        (2 * comp["N2"], RelativisticCoulomb(7)),
        (2 * comp["O2"], RelativisticCoulomb(8)),
        (2 * comp["N2"], SeltzerBerger(7)),
        (2 * comp["O2"], SeltzerBerger(8)),
        
        [(comp["N2"], orb) for orb in ORBITALS["N2"]]...,
        [(comp["O2"], orb) for orb in ORBITALS["O2"]]...,
    
        # (2 * comp["N2"], Moller(7, tcut)),
        # (2 * comp["O2"], Moller(8, tcut))
    ]
        
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

end

# if abspath(PROGRAM_FILE) == @__FILE__
#     Swarm.main()
# end

