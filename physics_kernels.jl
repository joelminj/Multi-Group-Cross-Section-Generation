# Physics calculations: Doppler broadening, flux spectra

# Physical Constants for kernel
const k_B = 8.617333262145e-5    # Boltzmann constant in eV/K
const A_Pu239 = 239.0            # Atomic mass of Pu-239

"""
    weighting_spectrum(E, T_n=573.15; isotope="Pu239")
Computes a continuous 3-region flux spectrum:
1. Thermal region: Maxwell-Boltzmann distribution
2. Epithermal region: 1/E (Fermi) spectrum
3. Fast region: Watt fission spectrum (isotope-dependent)

Arguments:
- E: Incident neutron energy (eV)
- T_n: Neutron temperature (K), default 573.15 K (300°C)
- isotope: Fissioning isotope ("U235" or "Pu239")
"""
function weighting_spectrum(E, T_n=573.15; isotope="Pu239")
    E_th = 0.2       # Transition to 1/E region (eV)
    E_fast = 1e5     # Transition to Watt Fission region (eV)
    
    # 1. Thermal: Maxwell-Boltzmann
    if E < E_th
        return (E / (k_B * T_n)^2) * exp(-E / (k_B * T_n))
        
    # 2. Epithermal: 1/E (Fermi)
    elseif E < E_fast
        # Calculate continuity constant C
        C = (E_th^2 / (k_B * T_n)^2) * exp(-E_th / (k_B * T_n))
        return C / E
        
    # 3. Fast: Watt Fission Spectrum
    else
        # Watt spectrum parameters based on ENDF/B evaluations
        if isotope == "U235"
            a = 0.988  # MeV^-1
            b = 2.249  # MeV^-1
        elseif isotope == "Pu239"
            a = 0.966  # MeV^-1
            b = 2.842  # MeV^-1
        else
            error("Watt parameters not defined for isotope: $isotope")
        end
        
        E_MeV = E / 1e6
        watt = exp(-a * E_MeV) * sinh(sqrt(b * E_MeV))
        
        # Calculate continuity constant D
        C = (E_th^2 / (k_B * T_n)^2) * exp(-E_th / (k_B * T_n))
        phi_epi_fast = C / E_fast
        watt_fast = exp(-a * (E_fast/1e6)) * sinh(sqrt(b * (E_fast/1e6)))
        D = phi_epi_fast / watt_fast
        
        return D * watt
    end
end
