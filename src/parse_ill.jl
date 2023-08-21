"""
parse_ill.jl

Load and parse ILL format scan file to DataFrame for post-processing
"""


function read_ill(fs::String)
    numor::String = ""
    instr::String = ""
    column_start::Int = -1
    data_start::Int = -1
    open(fs) do file
        for (index, line) in enumerate(eachline(file))
            # if occursin("INSTR: ", line)
            #     instr = line[7:end]
            # end
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
    display(numor[2:end])
    DataFrame(CSV.File(fs, header=column_start, skipto=data_start,
                       delim=" ", ignorerepeated=true, stripwhitespace=true,
                       drop=[:PNT], types=Float64))
end

# using BenchmarkTools
# @benchmark read_ill("./eiger2022n003415.scn")
# fs = read_ill("./eiger2022n003415.scn")