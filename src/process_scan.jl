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
                               :INSTR => first,
                               renamecols=false), :bin_labels)
    grouped_df[!, :I] .= grouped_df[!, :CNTS] ./ grouped_df[!, ncol]
    grouped_df[!, :I_ERR] .= sqrt.(grouped_df[!, :CNTS]) ./ grouped_df[!, ncol]
    grouped_df
end


function add_scans(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                  ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                  binsize::Float64=0.005)::DataFrame
    df_all = vcat([parse_numor(data_prefix, numor=n, ncol=ncol) for n in numors]...,
                 cols=:union)
    added_numors = join(unique(df_all.NUMOR), "_")
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    df_add = combine(groupby(df_all, :bin_labels),
                    xcol => mean,
                    ycol => sum,
                    ncol => sum,
                    :INSTR => first,
                    renamecols=false)
    df_add = unique(df_add, :bin_labels)
    df_add[!, :NUMOR] .= added_numors
    df_add[!, :I] .= df_add[!, :CNTS] ./ df_add[!, ncol]
    df_add[!, :I_ERR] .= sqrt.(df_add[!, :CNTS]) ./ df_add[!, ncol]
    df_add
end


function sub_scans(df_bg::DataFrame, df_fg::DataFrame, xcol::Symbol;
                  ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                  binsize::Float64=0.005)
    bg, fg = copy(df_bg), copy(df_fg)
    minx, maxx = extrema(bg[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    for df in [bg, fg]
        df.bin_labels = cut(df[!, xcol], linear_bins, extend=false)
        df = combine(groupby(df, :bin_labels),
                   xcol => mean,
                   ycol => sum,
                   ncol => sum,
                   :NUMOR => .*,
                   :INSTR => first,
                   renamecols=false)
        df = unique(df, :bin_labels)
        df[!, :I] .= df[!, :CNTS] ./ df[!, ncol]
        df[!, :I_ERR] .= sqrt.(df[!, :CNTS]) ./ df[!, ncol]
    end
    gbg = groupby(bg, :bin_labels)
    gfg = groupby(fg, :bin_labels)
    # for df in gbg
    #     display(df)
    # end
    for df_bin in gbg
        try
            I_fg = fg[fg.bin_labels .== df_bin.bin_labels, :].I
            I_ERR_fg = fg[fg.bin_labels .== df_bin.bin_labels, :].I_ERR
        catch ex
            if isa(ex, LoadError)
                I_fg = 0.0
                I_fg = 0.0
            end
        end
        # I_bg = df_bin.I
        # I_ERR_bg = df_bin.I_ERR
        # display(I_fg)
    end
    bg, fg
end