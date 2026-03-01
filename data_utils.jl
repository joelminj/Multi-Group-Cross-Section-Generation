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
