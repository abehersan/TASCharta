"""
parse_ill.jl

Parse ILL format file to DataFrame for post-processing
"""


"""
    parse_ill_file(filepath::String; ncol::Symbol=:M1)::DataFrame

Parse raw ASCII file that contains an ILL-formatted scan.
`filepath` is a string that contains the relative path of the scan.
`ncol` is the normalization (monitor) column, passed as a `Symbol`.

Scan columns are read and returned in `DataFrame` format.
"""
function parse_ill_file(filepath::String; ncol::Symbol=:M1)::DataFrame
    numor::String = ""
    instr::String = ""
    column_start::Int = -1
    data_start::Int = -1
    open(filepath) do file
        for (index, line) in enumerate(eachline(file))
            if occursin("INSTR: ", line)
                instr *= line[7:end]
            end
            if occursin("FILE_: ", line)
                numor *= line[7:end]
            end
            if occursin("DATA_:", line)
                column_start = index + 1
                data_start = index + 2
                break
            end
        end
    end
    df::DataFrame = DataFrame(CSV.File(filepath, header=column_start, skipto=data_start,
                       delim=" ", ignorerepeated=true, stripwhitespace=true,
                       drop=[:PNT], types=Float64))
    df[!, :NUMOR] .= numor[2:end]
    df[!, :INSTR] .= instr[2:end]
    df[!, :I] .= df[!, :CNTS] ./ df[!, ncol]
    df[!, :I_ERR] .= sqrt.(df[!, :CNTS]) ./ df[!, ncol]
    df
end


"""
    parse_numor_ill(data_prefix::String; numor::Int64, ncol::Symbol=:M1)::DataFrame

Parse numor given a data prefix string.
`data_prefix` is a formattable string that contains the relative path to the
data scans. The formattable part of the string is the numor itself.
Usually of the form 'tasp2023n%06d.dat'.
`numor` is an integer to be replaced in the prefix.
`ncol` is the normalization (monitor) column.

Scan columns are read and returned in `DataFrame` format.
"""
function parse_numor_ill(data_prefix::String; numor::Int64, ncol::Symbol=:M1)::DataFrame
    parse_ill_file(Printf.format(Printf.Format(data_prefix), numor), ncol=ncol)
end


save_scan(savepath::String, df::DataFrame) = CSV.write(savepath, df, delim="\t")
