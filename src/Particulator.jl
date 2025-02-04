module Particulator
export LogLinRange, MultiPopulation, RelativisticCoulomb, totalcs, CollisionTable,
    NullCollision, ElectronState, PhotonState, Population, ORBITALS, advance!, remove_particle!,
    repack!, energy, nparticles, AbstractCollisionTracker, track, SeltzerBerger

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
include("particledefs.jl")
include("population.jl")
include("mixed_population.jl")
include("collisions.jl")
include("electron.jl")
include("photon.jl")
include("relativistic_coulomb.jl")
include("rbeb.jl")
include("seltzer.jl")
include("slow-electron.jl")
include("zhelezniak_photon.jl")
include("lxcat.jl")

end
