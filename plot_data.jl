
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
