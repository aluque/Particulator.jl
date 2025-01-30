module Particulator
export LogLinRange, MultiPopulation, RelativisticCoulomb, totalcs, CollisionTable,
    NullCollision, ElectronState, Population, ORBITALS, advance!, remove_particle!,
    repack!, energy, nparticles, AbstractCollisionTracker, track

using Base.Threads: Atomic, @threads, atomic_add!
using LinearAlgebra
using StaticArrays
using Interpolations
using Polyester
using StructArrays
using DocStringExtensions
using Distributions
using LambertW

import JSON

include("constants.jl")
const co = Constants

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
include("relativistic_coulomb.jl")
include("rbeb.jl")
include("slow-electron.jl")
include("zhelezniak_photon.jl")
include("lxcat.jl")

end
