function combination_typeaware(x::AbstractArray)
    if eltype(x) == String
        return first(x)
    elseif eltype(x) <: Real
        return mean(x)
    else
        return first(x)
    end
end


@doc raw"""
    normalize_counts!(df::DataFrame; ycol::Symbol, ncol::Symbol)::Nothing

Normalize the recorded counts per bin to the standard monitor.
`ycol` is the raw counts column and `ncol` the column of the normalization.

Columns `I` and `I_ERR` are written to `df` assuming Poisson counting statistics.
"""
function normalize_counts!(df::DataFrame; ycol::Symbol, ncol::Symbol)::Nothing
    df[!, :I] .= df[!, ycol] ./ df[!, ncol]
    df[!, :I_ERR] .= sqrt.(df[!, ycol]) ./ df[!, ncol]
    return nothing
end


@doc raw"""
    rebin_scan(df::DataFrame, xcol::Symbol, ycol::Symbol, binsize::Float64)::DataFrame

Rebin scan along the `xcol` column. The size of the bins is specified as a single
float `binsize`.
"""
function rebin_scan(df::DataFrame, xcol::Symbol, ycol::Symbol, dycol::Symbol, binsize::Float64)::DataFrame
    minx, maxx = extrema(df[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    df[!, :bin_labels] .= cut(df[!, xcol], linear_bins, extend=true)
    df_bin = unique(combine(groupby(df, :bin_labels),
                            xcol => mean,
                            ycol => mean,
                            dycol => x->sqrt.(sum(x .^2)) ./ length(x),
                            renamecols=false), :bin_labels)

    all_numors = join(unique(df[!, :NUMOR]), "_")
    df_bin[!, :NUMOR] .= all_numors
    return df_bin
end


@doc raw"""
    add_numors(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                    ycol::Symbol=:CNTS, ncol::Symbol=:M1, parse_func::Function,
                    bins::Float64=0.005)::DataFrame

Add numors specified in the `numors` vector of integers.
`data_prefix` is a formattable string that specifies the location and numbering
of the datafiles.
`parse_func` is a parsing function for each individual datafile.

Scans are first binned, added and the normalized intensity is calculated
from the `ycol` and `ncol` symbols.
"""
function add_numors(data_prefix::String, numors::Vector{Int64}, xcol::Symbol;
                    ycol::Symbol=:CNTS, ncol::Symbol=:M1, parse_func::Function,
                    bins::Float64=0.005)::DataFrame
    df_all = vcat([parse_func(data_prefix, numor=n) for n in numors]..., cols=:union)
    added_numors = join(unique(df_all.NUMOR), "_")
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = collect((minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins))
    # linear_bins = minx:bins:maxx
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    col_symbs = setdiff(propertynames(df_all), [xcol, ycol, ncol, :bin_labels])
    df_add = unique(combine(groupby(df_all, :bin_labels),
                                       col_symbs .=> combination_typeaware,
                                       xcol => mean,
                                       ycol => sum,
                                       ncol => sum,
                                       renamecols=false), :bin_labels)
    df_add[!, :NUMOR] .= added_numors
    normalize_counts!(df_add, ycol=ycol, ncol=ncol)
    return df_add
end


@doc raw"""
    add_scans(dfs::Vector{DataFrame}, xcol::Symbol;
                    ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                    bins::Float64=0.005)::DataFrame

Add scans ad given in the vector of `DataFrame` `dfs`.

Scans are first binned, added and the normalized intensity is calculated
from the `ycol` and `ncol` symbols.
"""
function add_scans(dfs::Vector{DataFrame}, xcol::Symbol;
                    ycol::Symbol=:CNTS, ncol::Symbol=:M1,
                    bins::Float64=0.005)::DataFrame
    df_all = vcat(dfs..., cols=:union)
    added_numors = join(unique(df_all.NUMOR), "_")
    minx, maxx = extrema(df_all[!, xcol])
    linear_bins = collect((minx-(minx%bins)-bins):bins:(maxx-(maxx%bins)+2*bins))
    # linear_bins = minx:bins:maxx
    df_all.bin_labels = cut(df_all[!, xcol], linear_bins, extend=true)
    col_symbs = setdiff(propertynames(df_all), [xcol, ycol, ncol, :bin_labels])
    df_add = unique(combine(groupby(df_all, :bin_labels),
                                       col_symbs .=> combination_typeaware,
                                       xcol => mean,
                                       ycol => sum,
                                       ncol => sum,
                                       renamecols=false), :bin_labels)
    df_add[!, :NUMOR] .= added_numors
    normalize_counts!(df_add, ycol=ycol, ncol=ncol)
    return df_add
end


@doc raw"""
    sub_scans(df_bg::DataFrame, df_fg::DataFrame, xcol::Symbol, ycol::Symbol, dycol::Symbol, binsize::Float64=0.005)::DataFrame

Subtract `df_bg` from `df_fg`.
Scans are first individually binned to common bins with size `binsize`.
Normalized intensities are then directly subtracted and the error is
correctly propagated in quadrature.
"""
function sub_scans(df_bg::DataFrame, df_fg::DataFrame,
                    xcol::Symbol, ycol::Symbol, dycol::Symbol,
                    binsize::Float64=0.005)::DataFrame
    bg = rebin_scan(df_bg, xcol, ycol, dycol, binsize)
    fg = rebin_scan(df_fg, xcol, ycol, dycol, binsize)
    df_sub = DataFrame([T[] for T in eltype.(eachcol(fg))], names(fg))
    for fg_pnt in eachrow(fg)
        bin = fg_pnt.bin_labels
        bg_pnt = filter(:bin_labels => ==(bin), bg)
        if isempty(bg_pnt)
            continue
        end
        bg_pnt = first(bg_pnt)
        sub_pnt = fg_pnt
        sub_pnt.I -= bg_pnt.I
        sub_pnt.I_ERR = sqrt(bg_pnt.I_ERR^2 + fg_pnt.I_ERR^2)
        sub_pnt.NUMOR = fg_pnt.NUMOR*"-"*bg_pnt.NUMOR
        push!(df_sub, sub_pnt)
    end
    return df_sub
end


@doc raw"""
    interp_datagrid(df::DataFrame; x::Symbol, y::Symbol, z::Symbol, dy::Float64=0.025, SDIG::Int=3, ATOL::Float64=1e-3)

Interpolate the contents of a `DataFrame` for heatmap plotting.
The data is only interpolated along one direction, namely `y` with Interpolation
step `dy`.
`SDIG` controls the digits for rounding of floats in the `x` and `y` axes.
`ATOL` controls the absolute tolerance for float comparison for indexing.
"""
function interp_datagrid(df::DataFrame; x::Symbol, y::Symbol, z::Symbol, dy::Float64=0.025, SDIG::Int=3, ATOL::Float64=1e-3)
    XX = sort(unique(round.(df[!, x], digits=SDIG)))
    YY = sort(unique(round.(df[!, y], digits=SDIG)))
    minY, maxY = extrema(YY)
    YYmod = minY:dy:maxY
    ZZ = Matrix{Float64}(undef, length(YYmod), length(XX))
    @inbounds for i in eachindex(XX)
        dfX = filter(r->isapprox(r[x], XX[i], atol=ATOL), df)
        sort!(dfX, y)
        ZZinterpolation = LinearInterpolation(dfX[!, y], dfX[!, z], extrapolation_bc=Line())
        ZZmod = ZZinterpolation(YYmod)
        ZZ[:, i] = ZZmod
    end
    return [XX, YYmod, ZZ]
end