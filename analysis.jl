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
