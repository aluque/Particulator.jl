module Particulator
export LogLinRange, MultiPopulation, RelativisticCoulomb, totalcs, CollisionTable,
    Electron, Positron, Photon,
    NullCollision, ElectronState, PositronState, PhotonState, Population, ORBITALS,
    nactives, advance!, remove_particle!, eachactive, init!, ChebyshevCollisionTable,
    repack!, kinenergy, nparticles,
    StateChangeOutcome, ReplaceParticleOutcome, ReplaceParticlePairOutcome, NewParticleOutcome,
    StateChangeOutcome, RemoveParticleOutcome, NullOutcome,
    RBEB, SeltzerBerger, PhotoElectric,
    BetheHeitler, PositronAnihilation, Compton, KleinNishinaCompton, Moller, Bhaba, BinaryIntervals,
    chebfit, speed,
    AbstractCallback, CombinedCallback, CollisionCounter, WallCallback, VoidCallback, RouletteCallback,
    ParticleCountCallback, PopulationTargetCallback, collision_table_from_processes,
    HomogeneousField, DoubleLayerField, ConfinedDoubleLayerField, run!, roulette!,
    ElectromagneticField, ContinuumLoss, ChebContinuumLoss, RK2Pusher, RestrictedPusher, NullPusher,
    NullForcing, CombinedForcing,
    RestrictedForcing


using Base.Threads: Atomic, @threads, atomic_add!
using StyledStrings
using LinearAlgebra
using StaticArrays
using Interpolations
using Polyester
using Accessors
using StructArrays
using DocStringExtensions
using Distributions
using DelimitedFiles
using LambertW
using BSplineKit
using Printf

import JSON


include("constants.jl")
const co = Constants

const DATA_DIR = joinpath(@__DIR__, "..", "data")

@template DEFAULT =
    """
    $(SIGNATURES)
    $(DOCSTRING)
    """

include("util.jl")
include("cheby.jl")
include("particledefs.jl")
include("callback.jl")
include("population.jl")
include("mixed_population.jl")
include("collision_table.jl")
include("collisions.jl")
include("electron.jl") 
include("positron.jl")
include("photon.jl")
include("pusher.jl")
include("relativistic_coulomb.jl")
include("rbeb.jl")
include("continuum.jl")
include("seltzer.jl")
include("static_sandia_data.jl")
include("atomic_shells.jl")
include("photo_electric.jl")
include("bethe_heitler.jl")
include("anihilation.jl")
include("compton.jl")
include("moller.jl")
include("bhaba.jl")
include("slow-electron.jl")
include("field.jl")
include("zhelezniak_photon.jl")
include("lxcat.jl")
include("run.jl")

end
