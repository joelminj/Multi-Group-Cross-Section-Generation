# RP Activity 2
# Author: Joel Minj

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

"""
    doppler_broaden_fast(E_target_arr, E_raw, XS_raw, T)
Optimized Doppler broadening using the SIGMA1 kernel.
Only integrates over ±4 Doppler widths around each target energy.

This is much faster than full broadening but still accurate for resonances.
"""
function doppler_broaden_fast(E_target_arr, E_raw, XS_raw, T)
    XS_broadened = zeros(length(E_target_arr))
    
    for (j, E) in enumerate(E_target_arr)
        Delta = sqrt(4 * k_B * T * E / A_Pu239)
        
        # Define integration window: ±8 Doppler widths
        dE_approx = 8 * Delta * sqrt(E)
        E_min = max(1e-5, E - dE_approx)
        E_max = E + dE_approx
        
        # Find array indices for this window
        i_min = searchsortedfirst(E_raw, E_min)
        i_max = searchsortedlast(E_raw, E_max)
        
        if i_max > i_min
            integral = 0.0
            # Trapezoidal integration within the window
            for k in 1:(i_max - i_min)
                Ep1, Ep2 = E_raw[i_min + k - 1], E_raw[i_min + k]
                XS1, XS2 = XS_raw[i_min + k - 1], XS_raw[i_min + k]
                
                term1_1 = exp(-((sqrt(E) - sqrt(Ep1))^2) / (Delta^2))
                term2_1 = exp(-((sqrt(E) + sqrt(Ep1))^2) / (Delta^2))
                y1 = XS1 * (Ep1 / E) * (term1_1 - term2_1)
                
                term1_2 = exp(-((sqrt(E) - sqrt(Ep2))^2) / (Delta^2))
                term2_2 = exp(-((sqrt(E) + sqrt(Ep2))^2) / (Delta^2))
                y2 = XS2 * (Ep2 / E) * (term1_2 - term2_2)
                
                integral += 0.5 * (Ep2 - Ep1) * (y2 + y1)
            end
            XS_broadened[j] = (1 / (Delta * sqrt(pi))) * integral
        else
            # Fallback: return unbroadened value via linear interpolation
            XS_broadened[j] = interp_linear(E, E_raw, XS_raw)
        end
    end
    
    return XS_broadened
end

