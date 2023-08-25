"""
ybbr3_elastic.jl

Load, bin, add and subtract YbBr3 elastic scans at TASP
"""


using TASCharta
using Plots


data_prefix = "../data/tasp2023n%06d.dat"
numors_200mK = [1375, 1376, 1377]
numors_850mK = [1384, 1385, 1386] # T=850mK
numors_10K = [1393, 1394, 1395]


"""
test to see if binning of a single scan (dframe) works as intended
    XXX: works!!!
"""
test = true
if test
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    tas1 = parse_numor(data_prefix, numor=numors_10K[1], ncol=:M1)
    plot!(fig, tas1.K, tas1.I, yerr=tas1.I_ERR, label="Original")
    display(tas1)

    minx, maxx = extrema(tas1[!, :K])
    bins = 0.015
    linear_bins = (minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins)
    # bin_scan!(tas1, :K, linear_bins)
    tas1b = bin_scan(tas1, :K, bins)
    scatter!(fig, tas1b.K, tas1b.I, yerr=tas1b.I_ERR, label="Binned")
    display(fig)
    display(tas1b)
end


"""
test of added scans and subsequent rebinning of the data
    XXX: works!!!
"""
test = false
if test
    tas_ub = add_scans(data_prefix, numors_10K, :K, bins=0.0051)
    tas_bb = add_scans(data_prefix, numors_10K, :K, bins=0.015)
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    plot!(fig, tas_ub.K, tas_ub.I, yerr=tas_ub.I_ERR, label="Unbinned added")
    scatter!(fig, tas_bb.K, tas_bb.I, yerr=tas_bb.I_ERR, label="Binned added")
    display(fig)
end


"""
test to see if subtracting two dataframes works as intended
    XXX: works!!!
"""
test = false
if test
    bg = add_scans(data_prefix, numors_10K, :K, bins=0.0051)
    fg = add_scans(data_prefix, numors_850mK, :K, bins=0.0051)
    sub = sub_scans(bg, fg, :K, bins=0.0051);
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    scatter!(fig, bg.K, bg.I, yerr=bg.I_ERR, label="BG")
    scatter!(fig, fg.K, fg.I, yerr=fg.I_ERR, label="FG")
    scatter!(fig, sub.K,sub.I, yerr=sub.I_ERR, label="SUB")
    display(fig)
end