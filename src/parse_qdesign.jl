function parse_qdesign_file(filepath::String)::DataFrame
    column_start = -1
    data_start = -1
    open(filepath) do file
        for (index, line) in enumerate(eachline(file))
            if occursin("[Data]", line)
                column_start = index + 1
                data_start = index + 2
                break
            end
        end
    end
    df = DataFrame(CSV.File(filepath, header=column_start, skipto=data_start,
                    delim=",", silencewarnings=true))
    return df
end