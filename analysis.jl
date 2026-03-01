# RP Activity 2
# Autor: Joel Minj


# Multigroup collapsing, resonance self-shielding, Westcott g-factor

"""
    perform_multigroup_analysis(data, E_tot_raw, XS_tot_raw; T_mod=573.15, sigma_0=50.0)

Performs 12-group self-shielded multigroup collapsing on all reaction channels.
Returns both infinite dilution and self-shielded cross sections, plus flux depression data.

Arguments:
- data: Dict of (name => (E, XS) tuples)
- E_tot_raw, XS_tot_raw: Total cross section arrays (for flux depression)
- T_mod: Moderator temperature (K)
- sigma_0: Background cross section (barns) for Narrow Resonance approximation

Returns: (mg_results_inf, mg_results_shielded, flux_profile_data)
"""
function perform_multigroup_analysis(data, E_tot_raw, XS_tot_raw; T_mod=573.15, sigma_0=50.0)
    # Physical Constants
    k_B = 8.617333262145e-5
    E_0 = 0.0253  # Thermal reference energy in eV
    T_0 = 293.6   # Reference temperature (20 C) in K
    
    # 12-Group Structure
    group_edges = [20e6, 1e6, 1e5, 1e4, 1e3, 100.0, 10.0, 1.0, 0.625, 0.3, 0.1, 0.025, 1e-5]
    num_groups = length(group_edges) - 1
    
    # Storage for results
    mg_results_inf = Dict{String, Vector{Float64}}()
    mg_results_shielded = Dict{String, Vector{Float64}}()
    flux_profile_E = Float64[]
    flux_profile_inf = Float64[]
    flux_profile_shielded = Float64[]
    
    # Process each reaction channel
    for (name, (E_raw, XS_raw)) in data
        println("\n==================================================")
        println("RESULTS FOR Pu-239: $name")
        println("==================================================")
        
        xs_inf_arr = zeros(num_groups)
        xs_shielded_arr = zeros(num_groups)
        
        # Loop over groups
        for g in 1:num_groups
            E_upper = group_edges[g]
            E_lower = group_edges[g+1]
            
            idx = findall(x -> E_lower <= x <= E_upper, E_raw)
            if isempty(idx)
                println("Group $g: No data in range.")
                continue
            end
            
            E_g = E_raw[idx]
            XS_g_raw = XS_raw[idx]
            
            # Doppler broaden at target temperature
            if E_upper <= 10000.0
                XS_g = doppler_broaden_fast(E_g, E_raw, XS_raw, T_mod)
                XS_tot_g = doppler_broaden_fast(E_g, E_tot_raw, XS_tot_raw, T_mod)
            else
                XS_g = XS_g_raw
                XS_tot_g = interp_linear.(E_g, Ref(E_tot_raw), Ref(XS_tot_raw))
            end
            
            # Calculate flux profiles
            flux_macro = weighting_spectrum.(E_g, T_mod, isotope="U235")
            flux_shielded = flux_macro .* (sigma_0 ./ (XS_tot_g .+ sigma_0))
            
            # Save flux data for Fission around resonance
            if name == "Fission" && E_upper <= 1.0 && E_lower >= 0.1
                append!(flux_profile_E, E_g)
                append!(flux_profile_inf, flux_macro)
                append!(flux_profile_shielded, flux_shielded)
            end
            
            # Collapse group constants
            sigma_avg_inf = trapz(E_g, XS_g .* flux_macro) / trapz(E_g, flux_macro)
            sigma_avg_sh = trapz(E_g, XS_g .* flux_shielded) / trapz(E_g, flux_shielded)
            
            xs_inf_arr[g] = sigma_avg_inf
            xs_shielded_arr[g] = sigma_avg_sh
            
            @printf("Group %02d [%.1e to %.1e eV] | Inf: %8.2f b | Shielded: %8.2f b\n",
                    g, E_upper, E_lower, sigma_avg_inf, sigma_avg_sh)
        end
        
        # Westcott g-factor for entire thermal region
        idx_th = findall(x -> 1e-5 <= x <= 1.0, E_raw)
        E_th = E_raw[idx_th]
        XS_th_broad = doppler_broaden_fast(E_th, E_raw, XS_raw, T_mod)
        
        flux_maxwellian = pure_maxwellian.(E_th, T_mod)
        num_max = trapz(E_th, XS_th_broad .* flux_maxwellian)
        den_max = trapz(E_th, flux_maxwellian)
        sigma_avg_maxwellian = num_max / den_max
        
        sigma_E0 = interp_linear(E_0, E_raw, XS_raw)
        g_factor = (sigma_avg_maxwellian / sigma_E0) * (2 / sqrt(pi)) * sqrt(T_mod / T_0)
        
        println("  ------------------------------------------------")
        println("  -> σ at 0.0253 eV (0K):       ", round(sigma_E0, digits=2), " barns")
        println("  -> Pure Maxwellian Avg:       ", round(sigma_avg_maxwellian, digits=2), " barns")
        println("  -> Westcott g-factor (300°C): ", round(g_factor, digits=4))
        println("  ------------------------------------------------\n")
        
        mg_results_inf[name] = xs_inf_arr
        mg_results_shielded[name] = xs_shielded_arr
    end
    
    return mg_results_inf, mg_results_shielded, (flux_profile_E, flux_profile_inf, flux_profile_shielded)
end
