module Particulator
export LogLinRange, MultiPopulation, RelativisticCoulomb, totalcs, CollisionTable,
    Electron, Positron, Photon,
    NullCollision, ElectronState, PositronState, PhotonState, Population, ORBITALS,
    advance!, remove_particle!, eachactive, init!, ChebyshevCollisionTable,
    repack!, energy, nparticles, AbstractCollisionTracker, track, SeltzerBerger, PhotoElectric,
    BetheHeitler, PositronAnihilation, Compton, Bhaba, BinaryIntervals, fit, speed

using Base.Threads: Atomic, @threads, atomic_add!
using LinearAlgebra
using StaticArrays
using Interpolations
using Polyester
using StructArrays
using DocStringExtensions
using Distributions
using DelimitedFiles
using LambertW
using BSplineKit

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
include("population.jl")
include("mixed_population.jl")
include("collision_table.jl")
include("collisions.jl")
include("electron.jl")
include("positron.jl")
include("photon.jl")
include("relativistic_coulomb.jl")
include("rbeb.jl")
include("seltzer.jl")
include("static_sandia_data.jl")
include("atomic_shells.jl")
include("photo_electric.jl")
include("bethe_heitler.jl")
include("anihilation.jl")
include("compton.jl")
include("bhaba.jl")
include("slow-electron.jl")
include("zhelezniak_photon.jl")
include("lxcat.jl")

end
