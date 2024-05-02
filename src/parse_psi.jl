@doc raw"""
    parse_psi_file(filepath::String)::DataFrame

Parse raw ASCII file that contains a PSI-formatted scan.
`filepath` is a string that contains the relative path of the scan.

Scan columns are read and returned in `DataFrame` format.
"""
function parse_file_psi(filepath::String)::DataFrame
    titleline::Int = 1
    paramline::Int = 2
    stepmonline::Int = 3
    cntsline::Int = 4
    cntserrline::Int = -1

    file = open(filepath)
    lines = readlines(file)

    instr = split(lines[titleline], ",")[1]

    params = split(lines[paramline], ",")
    lmbda = parse(Float64, strip(split(params[1], "=")[end]))
    tt = parse(Float64, strip(split(params[2], "=")[end]))

    steps_and_mon = split(split(lines[stepmonline], ",")[1])
    ttmin = parse(Float64, steps_and_mon[1])
    ttstep = parse(Float64, steps_and_mon[2])
    ttmax = parse(Float64, steps_and_mon[3])
    meanmon = parse(Float64, steps_and_mon[end])

    ttrange = ttmin:ttstep:ttmax
    ttlen = length(ttrange)
    blocklength = fld(ttlen, 10)
    cntserrline = Int(blocklength + stepmonline) + 1

    CNTS = Int64[]
    for line in lines[cntsline:cntserrline - 1]
        counts = parse.(Int64, strip.(split(line), [only(".")]))
        push!(CNTS, counts...)
    end

    ERRS = Float64[]
    for line in lines[cntserrline:cntserrline + blocklength - 1]
        err = parse.(Float64, strip.(split(line), [only(".")]))
        push!(ERRS, err...)
    end

    @assert length(CNTS) == length(ERRS) == ttlen

    df = DataFrame(
        INSTR=instr,
        WVLEN=lmbda,
        TT=tt,
        MON=meanmon,
        TTHETA=ttrange,
        CNTS=CNTS, ERRS=ERRS
    )

    return df
end