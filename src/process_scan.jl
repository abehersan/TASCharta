"""
process_scan.jl

Add, bin and subtract dataframes resulting from read ILL scan files
"""


function combination_typeaware(x::AbstractArray)
    if eltype(x) == String
        first(x)
    elseif eltype(x) <: Real
        mean(x)
    else
        first(x)
    end
end


"""
    bin_scan(df::DataFrame, xcol::Symbol, bins::Float64; ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame

Bin the `xcol` column of a scan `DataFrame`.
The bins are determined linearly from the minimum and maximum of the `xcol`
of the scan.
`bins` is a `Float64` that determines the spacing between the bins.

Neutron counts that land in the same bin are added.
Normalization (monitor) is added as well.
The normalized intensities are returned in `DataFrame` format.
"""
function bin_scan(df::DataFrame, xcol::Symbol, bins::Float64;
                 ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame
    minx, maxx = extrema(df[!, xcol])
    linear_bins = (minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins)
    df.bin_labels = cut(df[!, xcol], linear_bins, extend=true)
    col_symbs = setdiff(propertynames(df), [xcol, ycol, ncol, :bin_labels])
    grouped_df::DataFrame = unique(combine(groupby(df, :bin_labels),
                                          col_symbs .=> combination_typeaware,
                                          xcol => mean,
                                          ycol => sum,
                                          ncol => sum,
                                          renamecols=false), :bin_labels)
    grouped_df[!, :I] .= grouped_df[!, :CNTS] ./ grouped_df[!, ncol]
    grouped_df[!, :I_ERR] .= sqrt.(grouped_df[!, :CNTS]) ./ grouped_df[!, ncol]
    grouped_df
end


"""
    bin_scan(df::DataFrame, xcol::Symbol, bins::StepRangeLen; ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame

Bin the `xcol` column of a scan `DataFrame`.
The bins are determined linearly from the minimum and maximum of the `xcol`
of the scan.
`bins` is a `StepRangeLen` that determines the breaks of the bins.

Neutron counts that land in the same bin are added.
Normalization (monitor) is added as well.
The normalized intensities are returned in `DataFrame` format.
"""
function bin_scan(df::DataFrame, xcol::Symbol, bins::StepRangeLen;
                 ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame
    df.bin_labels = cut(df[!, xcol], bins, extend=true)
    col_symbs = setdiff(propertynames(df), [xcol, ycol, ncol, :bin_labels])
    grouped_df::DataFrame = unique(combine(groupby(df, :bin_labels),
                                          col_symbs .=> combination_typeaware,
                                          xcol => mean,
                                          ycol => sum,
                                          ncol => sum,
                                          renamecols=false), :bin_labels)
    grouped_df[!, :I] .= grouped_df[!, :CNTS] ./ grouped_df[!, ncol]
    grouped_df[!, :I_ERR] .= sqrt.(grouped_df[!, :CNTS]) ./ grouped_df[!, ncol]
    grouped_df
end


"""
    add_scans(data_prefix::String, numors::Vector{Int64}, xcol::Symbol; ycol::Symbol=:CNTS, ncol::Symbol=:M1, bins::Float64=0.005)::DataFrame

Bin and add (combine) scans.
`data_prefix` is a formatted string for the location of the datafiles.
`numors` is a vector of numors.
`xcol` and `ycol` determine the columns that will be binned and added,
`bins` determines the size of the linearly-spaced bins.

The added scans are returned as a `DataFrame`.
"""
function add_scans(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                  ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                  bins::Float64=0.005)::DataFrame
    df_all::DataFrame = vcat([parse_numor(data_prefix, numor=n, ncol=ncol) for n in numors]..., cols=:union)
    added_numors::String = join(unique(df_all.NUMOR), "_")
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = (minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins)
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    col_symbs = setdiff(propertynames(df_all), [xcol, ycol, ncol, :bin_labels])
    df_add::DataFrame = unique(combine(groupby(df_all, :bin_labels),
                                       col_symbs .=> combination_typeaware,
                                       xcol => mean,
                                       ycol => sum,
                                       ncol => sum,
                                       renamecols=false), :bin_labels)
    df_add[!, :NUMOR] .= added_numors
    df_add[!, :I] .= df_add[!, :CNTS] ./ df_add[!, ncol]
    df_add[!, :I_ERR] .= sqrt.(df_add[!, :CNTS]) ./ df_add[!, ncol]
    df_add
end


"""
    sub_scans(df_bg::DataFrame, df_fg::DataFrame, xcol::Symbol; bins::Float64=0.005, ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame

Bin and subtract scans.
`df_bg` and `df_fg` are `DataFrame`s containing the scans of the background
and the foreground respectively.
`xcol` and `ycol` determine the columns that will be binned and subtracted,
`bins` determines the size of the linearly-spaced bins.

The subtracted scans are returned as a `DataFrame`.
"""
function sub_scans(df_bg::DataFrame, df_fg::DataFrame, xcol::Symbol;
                  bins::Float64=0.005, ycol::Symbol=:CNTS, ncol::Symbol=:M1)::DataFrame
    minx, maxx = extrema(df_fg[!, xcol])
    linear_bins = (minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins)
    bg::DataFrame = bin_scan(df_bg, xcol, linear_bins, ycol=ycol, ncol=ncol)
    fg::DataFrame = bin_scan(df_fg, xcol, linear_bins, ycol=ycol, ncol=ncol)
    df_sub::DataFrame = DataFrame([T[] for T in eltype.(eachcol(fg))], names(fg))
    for fg_pnt in eachrow(fg)
        bin = fg_pnt.bin_labels
        bg_pnt = first(filter(:bin_labels => ==(bin), bg))
        if isempty(bg_pnt)
            continue
        end
        sub_pnt = fg_pnt
        sub_pnt.I -= bg_pnt.I
        sub_pnt.I_ERR = sqrt(bg_pnt.I_ERR^2 + fg_pnt.I_ERR^2)
        sub_pnt.NUMOR = fg_pnt.NUMOR*"-"*bg_pnt.NUMOR
        push!(df_sub, sub_pnt)
    end
    df_sub
end