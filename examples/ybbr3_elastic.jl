"""
ybbr3_elastic.jl

Load, bin, add and subtract YbBr3 elastic scans at TASP
"""


using SINQCharta
using DataFrames
using StatsPlots


"""
test to see if binning of a single scan (dframe) works as intended
    XXX: works!!!
"""
function single_scan_bin()
    data_prefix = "../data/20210050/tasp2023n%06d.dat"
    numors_10K = [1393, 1394, 1395]

    p = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    tas1 = parse_numor_ill(data_prefix, numor=numors_10K[1])
    normalize_counts!(tas1, ycol=:CNTS, ncol=:M1)
    @df tas1 plot!(p, :K, :I, yerr=:I_ERR, label="Unbinned", ls=:dash, shape=:circle)

    # minx, maxx = extrema(tas1[!, :K])
    # bins = 0.015 # also an option to bin with just a step size
    # linear_bins = minx:bins:maxx
    tas1b = rebin_scan(tas1, :K, :I, :I_ERR, 0.01)
    @df tas1b scatter!(p, :K, :I, yerr=:I_ERR, label="Binned")
    display(names(tas1b))
    display(p)
    return nothing
end

single_scan_bin()



# """
# test of added scans and subsequent rebinning of the data
#     XXX: works!!!
# """
# function scan_addition()
#     data_prefix = "./data/20210050/tasp2023n%06d.dat"
#     numors_10K = [1393, 1394, 1395]

#     tas_ub = add_scans(data_prefix, numors_10K, :K, parse_func=parse_numor_ill, bins=0.0051)
#     tas_bb = add_scans(data_prefix, numors_10K, :K, parse_func=parse_numor_ill, bins=0.015)
#     normalize_counts!(tas_ub, ycol=:CNTS, ncol=:M1)
#     normalize_counts!(tas_bb, ycol=:CNTS, ncol=:M1)

#     p = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
#     @df tas_ub plot!(p, :K, :I, yerr=:I_ERR, label="Unbinned added")
#     scatter!(p, tas_bb.K, tas_bb.I, yerr=tas_bb.I_ERR, label="Binned added")
#     p
# end


# """
# test to see if subtracting two dataframes works as intended
#     XXX: works!!!
# """
# function scan_subtraction()
#     data_prefix = "./data/20210050/tasp2023n%06d.dat"
#     numors_10K = [1393, 1394, 1395]
#     numors_200mK = [1375, 1376, 1377]

#     bg = add_scans(data_prefix, numors_10K, :K, parse_func=parse_numor_ill, bins=0.0051)
#     fg = add_scans(data_prefix, numors_200mK, :K, parse_func=parse_numor_ill, bins=0.0051)
#     normalize_counts!(bg, ycol=:CNTS, ncol=:M1)
#     normalize_counts!(fg, ycol=:CNTS, ncol=:M1)

#     sub = sub_scans(bg, fg, :K, bins=0.0051);
#     display(first(sub.NUMOR))
#     # normalize_counts!(sub, ycol=:CNTS, ncol=:M1)

#     p = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
#     @df bg scatter!(p, :K, :I, yerr=:I_ERR, label="BG")
#     @df fg scatter!(p, :K, :I, yerr=:I_ERR, label="FG")
#     @df sub scatter!(p, :K, :I, yerr=:I_ERR, label="SUB")
#     p
# end


# function test_set()
#     ps = [single_scan_bin(), scan_addition(), scan_subtraction()]
#     plot(ps...)
# end


# test_set()
