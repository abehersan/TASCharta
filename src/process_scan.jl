"""
process_scan.jl

Add, bin and subtract dataframes resulting from read ILL scan files
"""

# TODO: always add CNTS, I, and ncol, propagate error in I and compare with
# poissonian error in added CNTS/added ncol

function bin_scan(df::DataFrame, xcol::Symbol;
                 ycol::Symbol=:I, ncol::Symbol=:M1,
                 binsize::Float64=0.005)::DataFrame
    minx, maxx = extrema(df[!, xcol])
    linear_bins = minx:binsize:maxx
    df.bin_labels = cut(df[!, xcol], linear_bins, extend=true)
    grouped_df = combine(groupby(df, :bin_labels),
                        xcol => mean,
                        ycol => sum,
                        ncol => sum,
                        :CNTS => sum,
                        :NUMOR => .*,
                        renamecols=false)
    # grouped_df[!, :I] .= grouped_df[!, :CNTS] ./ grouped_df[!, ncol]
    # grouped_df[!, :I_ERR] .= sqrt.(grouped_df[!, :CNTS]) ./ grouped_df[!, ncol]
    grouped_df
end


function add_scans(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                  ycol::Symbol=:I, ncol::Symbol=:M1,
                  binsize::Float64=0.005)::DataFrame
    df_all = vcat([parse_numor(data_prefix, numor=n, ncol=ncol)
                  for n in numors]..., cols=:union)
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = minx:binsize:maxx
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    df_add = combine(groupby(df_all, :bin_labels),
                    xcol => mean,
                    ycol => sum,
                    ncol => sum,
                    :CNTS => sum,
                    :NUMOR => .*,
                    renamecols=false)
    df_add = unique(df_add, :bin_labels)
    # df_add[!, :I] .= df_add[!, :CNTS] ./ df_add[!, ncol]
    # df_add[!, :I_ERR] .= sqrt.(df_add[!, :CNTS]) ./ df_add[!, ncol]
    df_add
end


function sub_scans(data_prefix::String, numors_bg::Vector{Int64},
                  numors_fg::Vector{Int64}, xcol::Symbol;
                  ycol::Symbol=:I, ncol::Symbol=:M1,
                  binsize::Float64=0.005)::DataFrame
    df_bg = add_scans(data_prefix, numors_bg, xcol, ycol=ycol, ncol=ncol, binsize=binsize)
    df_fg = add_scans(data_prefix, numors_fg, xcol, ycol=ycol, ncol=ncol, binsize=binsize)
    sub_df = combine(groupby(vcat(df_bg, df_fg, cols=:union), :bin_labels),
                    xcol => mean,
                    :I => -,
                    :I_ERR => x->sqrt(sum(x .^2)),
                    :NUMOR => .*,
                    renamecols=false)
    sub_df
end