# Data loading and numerical integration functions

using CSV, DataFrames

"""
    load_janis_csv(filepath)
Load nuclear data from JANIS-formatted CSV file.
Skips header rows and extracts Energy and CrossSection columns.
"""
function load_janis_csv(filepath)
    df = CSV.read(filepath, DataFrame, skipto=4, header=["Energy", "CrossSection"], 
                  select=[1, 2], delim=';', types=Float64)
    dropmissing!(df)
    return df.Energy, df.CrossSection
end

"""
    trapz(x, y)
Trapezoidal numerical integration for discrete arrays using the trapezoid rule.
"""
function trapz(x, y)
    integral = 0.0
    for i in 1:(length(x)-1)
        integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
    end
    return integral
end

"""
    interp_linear(E_target, E_arr, XS_arr)
Linear interpolation to find cross-section at a specific energy.
Returns boundary values if target is outside array range.
"""
function interp_linear(E_target, E_arr, XS_arr)
    idx = searchsortedfirst(E_arr, E_target)
    if idx == 1 return XS_arr[1] end
    if idx > length(E_arr) return XS_arr[end] end
    
    E1, E2 = E_arr[idx-1], E_arr[idx]
    XS1, XS2 = XS_arr[idx-1], XS_arr[idx]
    return XS1 + (XS2 - XS1) * (E_target - E1) / (E2 - E1)
end

"""
    logspace(start, stop, n)
Create a logarithmically-spaced array from 10^start to 10^stop with n points.
"""
function logspace(start, stop, n)
    return 10 .^ range(start, stop, length=n)
end

