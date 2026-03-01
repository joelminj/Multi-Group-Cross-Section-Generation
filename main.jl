# RP Activity 2
# Author: Joel Minj

using CSV, DataFrames, Plots, LaTeXStrings, Printf

gr()

# Components

include("data_utils.jl")
include("physics_kernels.jl")
include("analysis.jl")
include("visualization.jl")


# Parameters


const E_THERMAL_MAX = 1.0          # 1 eV
const E_FAST_MIN = 1e6             # 1 MeV
const E_FAST_MAX = 20e6            # 20 MeV
const T_MOD = 573.15               # Moderator Temperature: 300°C in K
const SIGMA_0 = 50.0               # Background cross section (barns)


# Data files to analyze

# Define files to analyze
files = Dict(
    "94-Pu-239 MT1 (n,total).csv" => "Total Cross Section",
    "94-Pu-239 MT2 (n,elastic).csv" => "Elastic Scattering",
    "94-Pu-239 MT18 (n,fission).csv" => "Fission",
    "94-Pu-239 MT102 (n,γ).csv" => "Radiative Capture"
)

dir = @__DIR__()

# Load all datasets
function load_all_data(files, dir)
    data = Dict()
    for (file, label) in files
        filepath = joinpath(dir, file)
        if !isfile(filepath)
            println("File not found: $file")
            continue
        end
        energy, cross_section = load_janis_csv(filepath)
        data[label] = (energy, cross_section)
    end
    return data
end

data = load_all_data(files, dir)

# Plot 0: Cross Sections with Thermal/Fast Regions


println("\nGenerating Cross-Section Plot")
p_cs = plot_cross_sections(data, E_THERMAL_MAX=E_THERMAL_MAX, E_FAST_MIN=E_FAST_MIN, E_FAST_MAX=E_FAST_MAX)
display(p_cs)
# Plot 0: Main cross-section overview
savefig(p_cs, joinpath(dir, "00_cross_sections.png"))



# Compute Multigroup Collapsing and Flux Depression

println("\n--- Calculating Multigroup Cross Sections at 300C ---")

E_tot_raw, XS_tot_raw = data["Total Cross Section"]

mg_results_inf, mg_results_shielded, flux_profile_data = perform_multigroup_analysis(
    data, E_tot_raw, XS_tot_raw,
    T_mod=T_MOD,
    sigma_0=SIGMA_0
)

# Plots

# Plot 1: Multigroup Collapse
println("\nGenerating Multigroup Collapse Plot...")
group_edges = [20e6, 1e6, 1e5, 1e4, 1e3, 100.0, 10.0, 1.0, 0.625, 0.3, 0.1, 0.025, 1e-5]
p_mg = plot_multigroup_collapse(data, mg_results_inf, mg_results_shielded, group_edges)
if p_mg !== nothing
    display(p_mg)
    savefig(p_mg, joinpath(dir, "01_multigroup_collapse.png"))
end

# Plot 2: Flux Depression
println("\nGenerating Flux Depression Plot...")
flux_profile_E, flux_profile_inf, flux_profile_shielded = flux_profile_data

E_tot_plot_idx = findall(x -> 0.1 <= x <= 1.0, E_tot_raw)
E_tot_plot = E_tot_raw[E_tot_plot_idx]
XS_tot_plot = doppler_broaden_fast(E_tot_plot, E_tot_raw, XS_tot_raw, T_MOD)

p_flux = plot_flux_depression((flux_profile_E, flux_profile_inf, flux_profile_shielded), E_tot_plot, XS_tot_plot)
display(p_flux)
savefig(p_flux, joinpath(dir, "02_flux_depression.png"))

println("\nFinished!")
