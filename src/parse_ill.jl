"""
parse_ill.jl

Parse ILL format file to DataFrame for post-processing
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


function parse_numor(data_prefix::String;
                    numor::Int64, ncol::Symbol=:M1)::DataFrame
    parse_ill_file(Printf.format(Printf.Format(data_prefix), numor), ncol=ncol)
end


function save_scan(savepath::String, df::DataFrame)
    CSV.write(savepath, df, delim="\t")
end