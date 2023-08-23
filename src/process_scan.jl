"""
process_scan.jl

Add, bin and subtract dataframes resulting from read ILL scan files
"""


function bin_scan(df::DataFrame, xcol::Symbol;
                 ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                 binsize::Float64=0.005)::DataFrame
    minx, maxx = extrema(df[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    df.bin_labels = cut(df[!, xcol], linear_bins, extend=true)
    grouped_df = unique(combine(groupby(df, :bin_labels),
                               xcol => mean,
                               ycol => sum,
                               ncol => sum,
                               :NUMOR => .*,
                               renamecols=false), :bin_labels)
    grouped_df[!, :I] .= grouped_df[!, :CNTS] ./ grouped_df[!, ncol]
    grouped_df[!, :I_ERR] .= sqrt.(grouped_df[!, :CNTS]) ./ grouped_df[!, ncol]
    grouped_df
end


function add_scans(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                  ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                  binsize::Float64=0.005)::DataFrame
    df_all = vcat([parse_numor(data_prefix, numor=n, ncol=ncol)
                  for n in numors]..., cols=:union)
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    df_add = combine(groupby(df_all, :bin_labels),
                    xcol => mean,
                    ycol => sum,
                    ncol => sum,
                    :NUMOR => .*,
                    renamecols=false)
    df_add = unique(df_add, :bin_labels)
    df_add[!, :I] .= df_add[!, :CNTS] ./ df_add[!, ncol]
    df_add[!, :I_ERR] .= sqrt.(df_add[!, :CNTS]) ./ df_add[!, ncol]
    df_add
end


function sub_scans(data_prefix::String, df_bg::DataFrame, df_fg::DataFrame,
                  xcol::Symbol; ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                  binsize::Float64=0.005)
    minx, maxx = extrema(df_bg[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    for df in [df_bg, df_fg]
        df.bin_labels = cut(df[!, xcol], linear_bins, extend=false)
        df= combine(groupby(df, :bin_labels),
                   xcol => mean,
                   ycol => sum,
                   ncol => sum,
                   :NUMOR => .*,
                   renamecols=false)
        df = unique(df, :bin_labels)
        df[!, :I] .= df[!, :CNTS] ./ df[!, ncol]
        df[!, :I_ERR] .= sqrt.(df[!, :CNTS]) ./ df[!, ncol]
    end
    df_bg, df_fg
end