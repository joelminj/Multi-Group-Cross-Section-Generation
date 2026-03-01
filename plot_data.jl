# RP Activity 2
# Author: Joel Minj

# Plotting functions for data


using Plots, LaTeXStrings

"""
    plot_cross_sections(data; E_THERMAL_MAX=1.0, E_FAST_MIN=1e6, E_FAST_MAX=20e6)

Creates a comprehensive log-log plot of all cross sections with thermal/fast regions marked.
"""
function plot_cross_sections(data; E_THERMAL_MAX=1.0, E_FAST_MIN=1e6, E_FAST_MAX=20e6)
    p = plot(
        xscale=:log10,
        yscale=:log10,
        xlabel="Incident Neutron Energy (eV)",
        ylabel="Cross Section (barns)",
        title="Pu-239 Nuclide Cross Sections",
        legend=:topright,
        grid=true,
        minorgrid=true,
        size=(900, 600),
        linewidth=1.5
    )
    
    for (name, (E, XS)) in data
        plot!(p, E, XS, label=name)
    end
    
    vspan!(p, [1e-5, E_THERMAL_MAX], color=:red, alpha=0.1, label="Thermal Region")
    vspan!(p, [E_FAST_MIN, E_FAST_MAX], color=:blue, alpha=0.1, label="Fast Region")
    
    return p
end

"""
    plot_multigroup_collapse(data, mg_results_inf, mg_results_shielded, group_edges)

Creates step plots comparing continuous data to 12-group collapsed cross sections.
Shows both infinite dilution and self-shielded multigroup results.
"""
function plot_multigroup_collapse(data, mg_results_inf, mg_results_shielded, group_edges)
    if !haskey(data, "Fission") || !haskey(mg_results_inf, "Fission") || !haskey(mg_results_shielded, "Fission")
        return nothing
    end
    
    p = plot(
        xscale=:log10, yscale=:log10,
        xlabel="Incident Neutron Energy (eV)",
        ylabel="Cross Section (barns)",
        title="Pu-239 Fission: Continuous vs Multigroup (Infinite Dilution & Self-Shielded)",
        legend=:topright,
        grid=true,
        minorgrid=true,
        size=(1000, 650),
        left_margin=8Plots.mm,
        bottom_margin=5Plots.mm
    )
    
    # Continuous data
    E_cont, XS_cont = data["Fission"]
    plot!(p, E_cont, XS_cont, label="Continuous (0K)", color=:lightgray, linewidth=1.5)
    
    # Infinite dilution multigroup step plot
    edges_asc = reverse(group_edges)
    xs_inf_asc = reverse(mg_results_inf["Fission"])
    xs_inf_step = push!(copy(xs_inf_asc), xs_inf_asc[end])
    
    plot!(p, edges_asc, xs_inf_step, st=:steppost,
          label="Infinite Dilution (300°C)", color=:red, linewidth=2.5)
    
    # Self-shielded multigroup step plot
    xs_sh_asc = reverse(mg_results_shielded["Fission"])
    xs_sh_step = push!(copy(xs_sh_asc), xs_sh_asc[end])
    
    plot!(p, edges_asc, xs_sh_step, st=:steppost,
          label="Self-Shielded σ₀=50b (300°C)", color=:blue, linewidth=2.5, linestyle=:dash)
    
    return p
end

"""
    plot_flux_depression(flux_profile_data, E_tot_plot, XS_tot_plot, T_mod=573.15)

Visualizes the Narrow Resonance flux depression effect at the 0.3 eV Pu-239 resonance.
"""
function plot_flux_depression(flux_profile_data, E_tot_plot, XS_tot_plot)
    flux_profile_E, flux_profile_inf, flux_profile_shielded = flux_profile_data
    
    # Sort by energy
    sort_idx = sortperm(flux_profile_E)
    E_plot = flux_profile_E[sort_idx]
    flux_inf_plot = flux_profile_inf[sort_idx]
    flux_sh_plot = flux_profile_shielded[sort_idx]
    
    p = plot(
        E_plot, flux_inf_plot,
        label="Macro Flux (Infinite Dilution)",
        color=:blue, linewidth=2.5,
        xlabel="Neutron Energy (eV)",
        ylabel="Relative Flux φ(E)",
        title="Neutron Flux Depression at the 0.3 eV Pu-239 Resonance",
        legend=:topleft, grid=true, size=(900, 600),
        xscale=:log10,
        left_margin=8Plots.mm,
        bottom_margin=5Plots.mm
    )
    
    plot!(p, E_plot, flux_sh_plot,
          label="Shielded Flux (σ₀ = 50 b)",
          color=:red, linewidth=2.5)
    
    # Add total cross section for reference
    scale_factor = maximum(flux_inf_plot) / maximum(XS_tot_plot)
    plot!(p, E_tot_plot, XS_tot_plot .* scale_factor,
          label="Total Cross Section (Scaled)",
          color=:black, linestyle=:dash, alpha=0.5, linewidth=2)
    
    return p
end

