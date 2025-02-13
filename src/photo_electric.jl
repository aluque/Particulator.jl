#= Photoelectric effect.
Implemented following GENAT4 and PENELOPE.
Total cross sections are obtained from the SANDIA data.
Angular sampling follows PENELOPE.
=#

struct PhotoElectric{T}
    Z::Int
    left_energy::Vector{T}
    coeffs::Matrix{T}


    # Note that we take a single binding energy for each element. We consider only K-shell ionization
    # and O2 and N2 have K-shell binding energy below 1 keV. For lower photon energies or for high-Z
    # elements we may have to select different binding energies  depending on the photon E but we
    # avoid that complication.
    binding::Vector{T}

    function PhotoElectric(T::Type, Z)
        s = (co.kilo * co.eV) .^ (1:4)

        A = Z / Z_TO_A_RATIO[Z]
        nintervs = NUMBER_OF_INTERVALS[Z]
        start = SANDIA_TABLE_INDEX[Z]
        data = SANDIA_TABLE[start:start+nintervs-1]
        
        left_energy = map(first, data) .* (co.kilo * co.eV)
        coeffs = hcat([(collect(d[2:end]) .* s) for d in data]...) .* (co.centi^2 * A / co.N_A) 
        binding = binding_energies(Z)
        
        return new{T}(Z, left_energy, coeffs, binding)
    end

    PhotoElectric(Z) = PhotoElectric(Float64, Z)
end


function collide(pe::PhotoElectric, photon::PhotonState{T}, eng) where T
    mc2 = co.electron_mc2
    electron_energy = sample_electron_energy(pe, eng)

    cosθ = sample_electron_cos_theta(pe, electron_energy)
    ϕ = 2π * rand()

    γ = 1 + electron_energy / mc2
    β = sqrt(1 - 1 / γ^2)
    v = turn(photon.p, cosθ, ϕ, co.c * β)

    ReplaceParticleOutcome(photon, ElectronState{T}(photon.x, v, photon.w, photon.t))
end

function totalcs(pe::PhotoElectric, eng)
    (;coeffs, left_energy) = pe
    j = searchsortedlast(left_energy, eng)    
    return sum(i -> (eng^-i) * coeffs[i, j], 1:4)
end

function sample_electron_energy(pe::PhotoElectric, photon_energy)
    (;binding) = pe
    
    local b
    for b1 in binding
        b = b1
        if photon_energy > b1
            break
        end
    end

    @assert photon_energy > b

    electron_energy = photon_energy - b

    return electron_energy
end

function sample_electron_cos_theta(::PhotoElectric, electron_energy)
    mc2 = co.electron_mc2

    # Sample outgoing angle with Sauter-Gavrila distribution.  We follow the PENELOPE manual (2015),
    # p. 55, eqs 2.6-2.11
    γ = 1 + electron_energy / mc2
    β = sqrt(1 - 1 / γ^2)
    A = 1 / β - 1
    g(ν) = (2 - ν) * (1 / (A + ν) + β * γ * (γ - 1) * (γ - 2) / 2)

    local ν
    accept = false
    while !accept
        ξ = rand()
        ν = 2A / ((A + 2)^2 - 4ξ) * (2ξ + (A + 2) * sqrt(ξ))
        
        ξ1 = rand()
        accept = ξ1 * g(0)  < g(ν)
    end
    
    return 1 - ν
end

