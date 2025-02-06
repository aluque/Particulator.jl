struct PositronAnihilation{T}
    # Z is allowed to be non-int to allow the use of effective Z, since everything is linear in Z
    Z::T
end

function totalcs(a::PositronAnihilation, eng)
    (;Z) = a
    mc2 = co.electron_mass * co.c^2

    γ = 1 + eng / mc2
    A = ((γ^2 + 4γ + 1) / (γ^2 - 1) * log(γ + sqrt(γ^2 - 1)) - (γ + 3) / sqrt(γ^2 - 1)) / (γ + 1)
    return Z * π * r_e^2 * A
end


"""
Sample photon energies as fraction of the initial positron TOTAL energy in a e+ + e-
anihilation event. The input `eng` is the positron KINETIC energy.
GEANT4 Phys. Ref. Manual v. 11.3 p. 146.
"""
function sample_secondary_energy_fraction(::PositronAnihilation, eng)
    mc2 = co.electron_mass * co.c^2

    γ = 1 + eng / mc2
    p = sqrt(T * (T + 2mc2)) / co.c

    ϵmax = (1 + sqrt((γ - 1) / (γ + 1)))
    ϵmin = (1 - sqrt((γ - 1) / (γ + 1)))
    # α = log(ϵmax / ϵmin)

    accept = false
    while !accept
        ϵ = ϵmin * (ϵmax / ϵmin)^rand()

        g = 1 - ϵ + (2γ * ϵ - 1) / (ϵ * (γ + 1)^2)
        accept = rand() < g
    end

    Etot = eng + 2mc2
    return ϵ
end

function secondary_cos_theta(::PositronAnihilation, eng, ϵ)
    γ = 1 + eng / mc2
    cosθ = (ϵ * (γ + 1) - 1) / (ϵ * sqrt(γ^2 - 1))

    return cosθ
end

