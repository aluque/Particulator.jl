struct Compton
    Z::Int
end

function collide(c::Compton, photon::PhotonState{T}, eng) where T
    mc2 = co.electron_mc2
    m2c2 = co.electron_mc^2

    E1, cosθ = sample_secondary_energy_and_cos_theta(c, eng)
    ϕ = 2π * rand()

    p1norm = E1 / co.c

    # Scattered photon momentum
    pg = turn(photon.p, cosθ, ϕ, p1norm)

    # Electron momentum from conservation
    pe = photon.p - pg

    pe2 = dot(pe, pe)
    β = sqrt(pe2 / (pe2 + m2c2))
    ve = co.c * β * pe / sqrt(pe2)
    
    photon1  = PhotonState{T}(photon.x, pg, photon.w, photon.t)
    electron = ElectronState{T}(photon.x, ve, photon.w, photon.t)

    return NewParticleOutcome(photon1, electron)
end


"""
Klein-Nishina cross-section following L&L 4, sect. 86.
"""
function klein_nishina(c::Compton, eng)
    (;Z) = c
    mc2 = co.electron_mc2

    # Classical electron radius, about 2.8e-15 m
    r_e = co.elementary_charge^2 / (co.electron_mass * co.c^2) / (4π * co.epsilon_0)

    x = 2 * eng / mc2

    sigma = (2π * r_e^2 * (1 / x) *
        ((1 - 4 / x - 8 / x^2) * log(1 + x) + 1 / 2 + 8 / x - 1 / (1 + x)^2 / 2))
    
    return Z * sigma
end


"""
Compute total cross-section of Compton scattering for energy `eng`.

Adapted from G4KleinNishinaCompton in GEANT4. 
For photon energies above ~100 keV this is the same as the Klein-Nishina cross section.
However it becomes much lower for lower energies. Note that only Z is used, neglecting information
about the atomic/molecular orbitals of electrons. However, this is likely much more accurate than the
Klein-Nishina cross-section.
"""
function totalcs(c::Compton, eng)
    (;Z) = c
    
    barn = 1e-28
    mc2 = co.electron_mc2
    
    a = 20.0
    b = 230.0
    c = 440.0

    d1 =  2.7965e-1 * barn
    d2 = -1.8300e-1 * barn
    d3 =  6.7527    * barn
    d4 = -1.9798e+1 * barn
    e1 =  1.9756e-5 * barn
    e2 = -1.0205e-2 * barn
    e3 = -7.3913e-2 * barn
    e4 =  2.7079e-2 * barn
    f1 = -3.9178e-7 * barn
    f2 =  6.8241e-5 * barn
    f3 =  6.0480e-5 * barn
    f4 =  3.0274e-4 * barn

    p1Z = Z * (d1 + e1 * Z + f1 * Z^2)
    p2Z = Z * (d2 + e2 * Z + f2 * Z^2)
    p3Z = Z * (d3 + e3 * Z + f3 * Z^2)
    p4Z = Z * (d4 + e4 * Z + f4 * Z^2)

    T0 = 15.0e3 * co.eV
    if Z <= 1
        T0 = 40.0 * keV
    end

    X = max(eng, T0) / mc2
    sigma = p1Z * log(1 + 2X) / X +
        (p2Z + p3Z * X + p4Z * X^2) / (1 + a * X + b * X^2 + c * X^3)

    # Modification for low energy (special case for Hydrogen)
    if eng < T0
        dT0 = 1e3 * co.eV
        X = (T0 + dT0) / mc2
        sigma1 = p1Z * log(1.0 + 2.0 * X) / X +
            (p2Z + p3Z * X + p4Z * X^2) / (1.0 + a * X + b * X^2 + c * X^3)
        c1 = -T0 * (sigma1 - sigma) / (sigma * dT0)
        c2 = 0.150
        if Z > 1.5
            c2 = 0.375 - 0.0556 * log(Z)
        end
        y = log(eng / T0)
        sigma *= exp(-y * (c1 + c2 * y))
    end

    return sigma
end

"""
Sample the energy of the secondaries and the scattering angle in a Compton collision.
Returns the energy of the scattered photon and the cosine of the angle of the scattered photon.
Based on the GEANT4 Phys. Ref. Manual v. 11.3 p. 41 and G4KleinNishinaCompton.cc.
"""
function sample_secondary_energy_and_cos_theta(::Compton, eng)
    mc2 = co.electron_mc2
    ϵ0 = mc2 / (mc2 + 2 * eng)
    α1 = -log(ϵ0)
    α2 = (1 - ϵ0^2) / 2

    local t, ϵ
    accept = false
    while !accept
        if rand() < α1 / (α1 + α2)
            ϵ = exp(-rand() * α1)
        else
            ϵ = sqrt(ϵ0^2 + (1 - ϵ0^2) * rand())
        end

        # t = 1 - cosθ
        t = mc2 * (1 - ϵ) / (ϵ * eng)
        g = (1 - ϵ / (1 + ϵ^2) * t * (2 - t))
        accept = rand() < g
    end    

    cosθ = 1 - t
    E1 = ϵ * eng

    return E1, cosθ
end
