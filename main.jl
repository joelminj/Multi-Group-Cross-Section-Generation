using CSV, DataFrames, Plots, LaTeXStrings, Printf

gr()


# Parameters


const E_THERMAL_MAX = 1.0          # 1 eV
const E_FAST_MIN = 1e6             # 1 MeV
const E_FAST_MAX = 20e6            # 20 MeV
const T_MOD = 573.15               # Moderator Temperature: 300°C in K
const SIGMA_0 = 50.0               # Background cross section (barns)
