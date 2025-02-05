struct BetheHeitler
    Z::Int
end

function collide(bh::BetheHeitler, photon::PhotonState{T}, eng) where T
    mc2 = co.electron_mass * co.c^2

    pkin, ekin = sample_secondary_energy(bh, eng)
    ϕ = 2π * rand()

    # Deal with the electron
    cosθ = sample_cos_theta(bh, ekin)
    γ = 1 + ekin / mc2
    β = (γ^2 - 1) / γ^2
    v_e = turn(photon.p, cosθ, ϕ, co.c * β)

    # Now the positron
    cosθ = sample_cos_theta(bh, pkin)
    γ = 1 + pkin / mc2
    β = (γ^2 - 1) / γ^2
    v_p = turn(photon.p, cosθ, ϕ, co.c * β)

    electron = ElectronState{T}(photon.x, v_e, photon.w, nextcoll(), photon.t, photon.active)
    positron = PositronState{T}(photon.x, v_p, photon.w, nextcoll(), photon.t, photon.active)
    
    return ReplaceParticlePairOutcome(photon, electron, positron)                                      
end

# Adapted from G4BetheHeitlerModel.cc
function totalcs(bh::BetheHeitler, eng)
    (;Z) = bh
    mc2 = co.electron_mass * co.c^2
    if eng < 2mc2
        # No pair production if gamma is low-energy
        return zero(eng)
    end

    microbarn = 1e-34

    gamma_energy_limit = 1.5e6 * co.eV;

    ## set coefficients a, b c
    a0 =  8.7842e+2 * microbarn
    a1 = -1.9625e+3 * microbarn 
    a2 =  1.2949e+3 * microbarn
    a3 = -2.0028e+2 * microbarn 
    a4 =  1.2575e+1 * microbarn 
    a5 = -2.8333e-1 * microbarn
    
    b0 = -1.0342e+1 * microbarn
    b1 =  1.7692e+1 * microbarn
    b2 = -8.2381    * microbarn
    b3 =  1.3063    * microbarn
    b4 = -9.0815e-2 * microbarn
    b5 =  2.3586e-3 * microbarn
    
    c0 = -4.5263e+2 * microbarn
    c1 =  1.1161e+3 * microbarn 
    c2 = -8.6749e+2 * microbarn
    c3 =  2.1773e+2 * microbarn 
    c4 = -2.0467e+1 * microbarn
    c5 =  6.5372e-1 * microbarn

    eng1 = max(eng, gamma_energy_limit)
    
    x = log(eng1 / mc2)
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x

    F1 = a0 + a1 * x + a2 * x2 + a3 * x3 + a4 * x4 + a5 * x5
    F2 = b0 + b1 * x + b2 * x2 + b3 * x3 + b4 * x4 + b5 * x5
    F3 = c0 + c1 * x + c2 * x2 + c3 * x3 + c4 * x4 + c5 * x5

    sigma = (Z + 1) * (F1 * Z + F2 * Z * Z + F3)
    
    if eng < gamma_energy_limit
        sigma = sigma * ((eng - 2mc2) / (gamma_energy_limit - 2mc2))^2
    end    

    return sigma
end

"""
Sample energy of secundary from primary (gamma) kinetic energy `eng` using the data in `bh`.
"""
function sample_secondary_energy(bh::BetheHeitler, eng)
    (;Z) = bh
    mc2 = co.electron_mass * co.c^2

    ϵ0 = mc2 / eng
    @assert ϵ0 < 0.5 "Gamma energy must be more than 2 * 511 keV to produce pairs"
    local ϵ
    
    if eng < 2e6 * co.eV
        ϵ = ϵ0 + (0.5 - ϵ0) * rand()
    else
        # Called deltaFactor in GEANT4.
        δ0 = 136 * ϵ0 / Z^(1/3)
        FZ = 8 * log(Z) / 3
        if eng > 50e6 * co.eV
            FZ += 8 * _fc(Z)
        end

        δmin = 4 * δ0
        δmax = exp((42.24 - FZ) / 8.368) - 0.952

        ϵp = (1 - sqrt(1 - δmin / δmax)) / 2
        ϵmin = max(ϵ0, ϵp)
        ϵrange = 0.5 - ϵmin
        
        F10, F20 = _screen_function12(δmin)
        F10 -= FZ
        F20 -= FZ

        NF1 = max(F10 * ϵrange^2, zero(F10))
        NF2 = max(1.5 * F20, zero(F20))
        NC = NF1 / (NF1 + NF2)
        
        accept = false
        while !accept
            if NC > rand()
                ϵ = 0.5 - ϵrange * rand()^(1/3)
                δ = δ0 / (ϵ * (1 - ϵ))
                accept = rand() < (_screen_function1(δ) - FZ) / F10
            else
                ϵ = ϵmin + ϵrange * rand()
                δ = δ0 / (ϵ * (1 - ϵ))
                accept = rand() < (_screen_function2(δ) - FZ) / F20
            end
        end
    end

    # Select e- / e+ randomly
    if rand(Bool)
        etotenergy = (1 - ϵ) * eng
        ptotenergy = ϵ * eng
    else
        ptotenergy = (1 - ϵ) * eng
        etotenergy = ϵ * eng
    end        

    # Find kinetic energies
    ekin = max(zero(etotenergy), etotenergy - mc2)
    pkin = max(zero(ptotenergy), ptotenergy - mc2)

    return (ekin, pkin)
end


sample_cos_theta(::BetheHeitler, T) = sample_modified_tsai_cos_theta(T)


""" Coulomb correction for > 50 MeV. """
function _fc(Z)
    αZ = co.fine_structure
    αZ2 = αZ * αZ
    αZ4 = αZ2 * αZ2
    αZ6 = αZ4 * αZ2

    f1 = 1 / (1 + αZ2) + 0.20206 − 0.0369 * αZ2 + 0.0083 * αZ4 − 0.0020 * αZ6
    return f1 * αZ2
end

function _screen_function12(delta)
    if delta > 1.4
        f1 = 42.038 - 8.29 * log(delta + 0.958)
        f2 = f1
    else
        f1 = 42.184 - delta * (7.444 - 1.623 * delta)
        f2 = 41.326 - delta * (5.848 - 0.902 * delta)
    end
    
    return (f1, f2)
end

function _screen_function1(delta)
    if delta > 1.4
        f1 = 42.038 - 8.29 * log(delta + 0.958)
    else
        f1 = 42.184 - delta * (7.444 - 1.623 * delta)
    end
    
    return f1
end

function _screen_function2(delta)
    if delta > 1.4
        f2 = 42.038 - 8.29 * log(delta + 0.958)
    else
        f2 = 41.326 - delta * (5.848 - 0.902 * delta)
    end

    return f2
end
