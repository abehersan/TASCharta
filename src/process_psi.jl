@doc raw"""
    rebin_psi(df::DataFrame, xcol::Symbol, ycol::Symbol, dycol::Symbol, binsize::Float64)::DataFrame

Rebin PSI-format data along the `xcol` column.
The size of the bins is specified as a single float `binsize`.
"""
function rebin_psi(df::DataFrame, xcol::Symbol, ycol::Symbol, dycol::Symbol, binsize::Float64)::DataFrame
    minx, maxx = extrema(df[!, xcol])
    linear_bins = (minx-(minx%binsize)-binsize):binsize:(maxx-(maxx%binsize)+2*binsize)
    df[!, :bin_labels] .= cut(df[!, xcol], linear_bins, extend=true)
    df_bin = unique(combine(groupby(df, :bin_labels),
                            xcol => mean,
                            ycol => mean,
                            dycol => x->sqrt.(sum(x .^2)) ./ length(x),
                            renamecols=false), :bin_labels)
    return df_bin
end


@doc raw"""
    add_psi(dfs::Vector{DataFrame}, xcol::Symbol,
            ycol::Symbol, dycol::Symbol, binsize::Float64)::DataFrame

Add all PSI dataframes given in the vector of `DataFrame`s `dfs`.

The data is first binned to common step in the `xcol`.
The resulting `DataFrame` is added and the uncertainty is propagated in quadrature.
"""
function add_psi(dfs::Vector{DataFrame}, xcol::Symbol, ycol::Symbol,
                dycol::Symbol, binsize::Float64)::DataFrame
    df_all = rebin_psi(vcat(dfs..., cols=:union), xcol, ycol, dycol, binsize)
    return df_all
end


@doc raw"""
    sub_psi(df_bg::DataFrame, df_fg::DataFrame, xcol::Symbol, ycol::Symbol, dycol::Symbol, binsize::Float64)::DataFrame

Calculate `df_fg` - `df_bg`.
The data are first individually binned to common bins with size `binsize`.
Counts are directly subtracted and the error is
correctly propagated in quadrature.
"""
function sub_psi(df_bg::DataFrame, df_fg::DataFrame,
                    xcol::Symbol, ycol::Symbol, dycol::Symbol,
                    binsize::Float64)::DataFrame
    bg = rebin_psi(df_bg, xcol, ycol, dycol, binsize)
    fg = rebin_psi(df_fg, xcol, ycol, dycol, binsize)
    df_sub = DataFrame([T[] for T in eltype.(eachcol(fg))], names(fg))
    for fg_pnt in eachrow(fg)
        bin = fg_pnt.bin_labels
        bg_pnt = filter(:bin_labels => ==(bin), bg)
        if isempty(bg_pnt)
            continue
        end
        bg_pnt = first(bg_pnt)
        sub_pnt = fg_pnt
        sub_pnt.CNTS -= bg_pnt.CNTS
        sub_pnt.ERRS = sqrt(bg_pnt.ERRS^2 + fg_pnt.ERRS^2)
        push!(df_sub, sub_pnt)
    end
    return df_sub
end