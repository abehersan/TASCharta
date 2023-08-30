"""
parse_zebra.jl

Parse ZEBRA point-detector ASCII files to DataFrame for post-processing
"""


"""
    parse_zebra_pointdet(filepath::String; ncol::Symbol=:Monitor1)::DataFrame

Parse raw ASCII file that contains a ZEBRA in point-detector mode.
`filepath` is a string that contains the relative path to the scan.
`ncol` is the normalization (monitor) column, passed as a `Symbol`.

Scan columns are read and returned in `DataFrame` format.
"""
function parse_zebra_pointdet(filepath::String; ncol::Symbol=:Monitor1)::DataFrame
    numor::String = ""
    instr::String = ""
    TS::Float64 = 0.0
    MF::Float64 = 0.0
    column_start::Int = -1
    data_start::Int = -1
    open(filepath) do file
        for (index, line) in enumerate(eachline(file))
            if occursin("instrument", line)
                instr *= split(line, "= ")[end]
            end
            if occursin("original_filename", line)
                numor *= replace(split(line, ".dat")[1][end-6:end], "0"=>"")
            end
            if occursin("temperature", line)
                TS = parse(Float64, split(line, "= ")[end])
            end
            if occursin("magnetic_field", line)
                MF = parse(Float64, split(line, "= ")[end])
            end
            if occursin("NP", line)
                column_start = index
                data_start = index + 1
                break
            end
        end
    end
    df::DataFrame = DataFrame(CSV.File(filepath, header=column_start, skipto=data_start,
                       delim=" ", ignorerepeated=true, stripwhitespace=true,
                       drop=[:NP], types=Float64, footerskip=1))
    cols = names(df)
    if in("B", cols)
        rename!(df, Dict(:B => :MF))
    else
        df[!, :MF] .= MF
    end
    if !in("TS", cols)
        df[!, :TS] .= TS
    end
    df[!, :NUMOR] .= numor[2:end]
    df[!, :INSTR] .= instr[1:end]
    df[!, :I] .= df[!, :Counts] ./ df[!, ncol]
    df[!, :I_ERR] .= sqrt.(df[!, :Counts]) ./ df[!, ncol]
    df
end


"""
    parse_numor_zebra(data_prefix::String; numor::Int64, ncol::Symbol=:M1)::DataFrame

Parse numor given a data prefix string.
`data_prefix` is a formattable string that contains the relative path to the
data scans. The formattable part of the string is the numor itself.
Usually of the form 'tasp2023n%06d.dat'.
`numor` is an integer to be replaced in the prefix.
`ncol` is the normalization (monitor) column.

Scan columns are read and returned in `DataFrame` format.
"""
function parse_numor_zebra(data_prefix::String; numor::Int64, ncol::Symbol=:Monitor1)::DataFrame
    parse_zebra_pointdet(Printf.format(Printf.Format(data_prefix), numor), ncol=ncol)
end
