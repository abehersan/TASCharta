"""
ybbr3_criticalfield.jl

Load, bin and add field-scans at ZEBRA to extract the critical field of YbBr3
"""


using SINQCharta
using StatsPlots
using DataFrames
default(size=(1080, 720))
default(margin=5Plots.mm)
default(linewidth=3.5)


function compare_ramps()
    data_prefix = "./data/20222632/zebra2023n%06d.dat"
    df_Bscan1 = add_scans(data_prefix, [5287, 5288, 5289], :MF, ycol=:Counts, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.1)
    df_Bscan2 = add_scans(data_prefix, [5299], :MF, ycol=:Counts, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.1)
    df_Bscan3 = add_scans(data_prefix, [5305, 5306], ycol=:Counts, :MF, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.1)
    df_Bscan_all = add_scans(data_prefix, [5287, 5288, 5289, 5299, 5305, 5306], ycol=:Counts, :MF, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.11)
    p = plot(xlabel="Magnetic Field [Tesla]", ylabel="Normalized Intensity",
            title="YbBr3, ZEBRA, T=65mK, Q=(1, 1, 0), B-scan", reuse=false)
    @df df_Bscan1 scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="NUMORS: $(first(:NUMOR))")
    @df df_Bscan2 scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="NUMORS: $(first(:NUMOR))")
    @df df_Bscan3 scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="NUMORS: $(first(:NUMOR))")
    @df df_Bscan_all scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="NUMORS: $(first(:NUMOR))")
    # display(df_Bscan1)
    ylims!(0, 5.0e-3)
    p
end


function field_ramps()
    data_prefix = "./data/20222632/zebra2023n%06d.dat"
    df_Bscan_110 = add_scans(data_prefix, [5287, 5288, 5289, 5299, 5305, 5306], ycol=:Counts, :MF, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.11)
    df_Bscan_100 = add_scans(data_prefix, [5307, 5308, 5309], ycol=:Counts, :MF, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.11)
    p = plot(xlabel="Magnetic Field [Tesla]", ylabel="Normalized Intensity",
            title="YbBr3, ZEBRA, T=65mK, B-scan", reuse=false)
    @df df_Bscan_110 scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="Q=(1, 1, 0), NUMORS: $(first(:NUMOR))")
    @df df_Bscan_100 scatter!(p, :MF, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="Q-(1, 0, 0), NUMORS: $(first(:NUMOR))")
    ylims!(0, 5.0e-3)
    p
end


function temperature_ramps()
    data_prefix = "./data/20222632/zebra2023n%06d.dat"
    df_Tscan = add_scans(data_prefix, [5310, 5311, 5312], ycol=:Counts, :T, ncol=:Monitor1, parse_func=parse_numor_zebra, bins=0.11)
    p = plot(xlabel="Temperature [K]", ylabel="Normalized Intensity",
            title="YbBr3, ZEBRA, B=11T, Q=(1, 1, 0), T-scan", reuse=false)
    @df df_Tscan scatter!(p, :T, :Counts ./ :Monitor1 , yerr=sqrt.(:Counts) ./ :Monitor1, label="Q=(1, 1, 0), B=11T")
    ylims!(0, 5.0e-3)
    p
end


function test_set()
    ps = [compare_ramps(), field_ramps(), temperature_ramps()]
    ps = [field_ramps(), temperature_ramps()]
    plot(ps...)
end

test_set()