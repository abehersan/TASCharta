"""
ybbr3_elastic.jl

Load, bin, add and subtract YbBr3 elastic scans at TASP
"""

using TASCharta
using Plots


data_prefix = "../data/tasp2023n%06d.dat"
numors_200mK = [1375, 1376, 1377]
numors_10K = [1393, 1394, 1395]


test = false
if test
    tas1 = parse_numor(data_prefix, numor=numors_10K[1], ncol=:M1)
    tas1b = bin_scan(tas1, :K, binsize=0.015)
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    scatter!(fig, tas1.K, tas1.I, yerr=tas1.I_ERR, label="Original")
    scatter!(fig, tas1b.K, tas1b.I, yerr=tas1b.I_ERR, label="Binned")
    display(fig)
    display(tas1)
    display(tas1b)
end


test = false
if test
    tas_ub = add_scans(data_prefix, numors_10K, :K, binsize=0.005)
    tas_bb = add_scans(data_prefix, numors_10K, :K, binsize=0.015)
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    scatter!(fig, tas_ub.K, tas_ub.I, yerr=tas_ub.I_ERR, label="Unbinned added")
    scatter!(fig, tas_bb.K, tas_bb.I, yerr=tas_bb.I_ERR, label="Binned added")
    display(fig)
end


test = true
if test
    bg = add_scans(data_prefix, numors_10K, :K, binsize=0.0051)
    fg = add_scans(data_prefix, numors_200mK, :K, binsize=0.0051)
    bg, fg = sub_scans(data_prefix, bg, fg, :K, binsize=0.015)
    fig = plot(xlabel="K [r.l.u.]", ylabel="Normalized Intensity")
    scatter!(fig, bg.K, bg.I, yerr=bg.I_ERR, label="BG")
    scatter!(fig, fg.K, fg.I, yerr=fg.I_ERR, label="FG")
    display(fig)
    display(bg)
    display(fg)
end