module Particulator

using Base.Threads: Atomic, @threads, atomic_add!
using LinearAlgebra
using StaticArrays
using Interpolations
using Polyester
using StructArrays
using DocStringExtensions
using Distributions

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
include("zphoton.jl")
include("lxcat.jl")

end
