module SINQCharta


using DataFrames
using Statistics
using CategoricalArrays
using CSV
using Printf

include("./parse_ill.jl")
export parse_file_ill
export parse_numor_ill

include("./parse_psi.jl")
export parse_file_psi

include("./process_psi.jl")
export rebin_psi
export add_psi
export sub_psi

include("parse_zebra.jl")
export parse_zebra_pointdet
export parse_numor_zebra

include("./process_ill.jl")
export normalize_counts!
export rebin_scan
export add_scans
export add_numors
export sub_scans

include("./parse_qdesign.jl")
export parse_qdesign_file

include("./process_qdesign.jl")
export get_chi
export get_mag
export get_hc
export make_massdict


end