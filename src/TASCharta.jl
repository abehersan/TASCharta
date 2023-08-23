module TASCharta

using DataFrames
using Statistics
using CategoricalArrays
using CSV
using Printf

include("parse_ill.jl")
export parse_ill_file
export parse_numor
export save_scan

include("process_scan.jl")
export bin_scan
export add_scans
export sub_scans

end
