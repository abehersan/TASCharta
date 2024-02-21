NA = 6.02214e23     # [1/mol]
muB = 9.274e-21     # [erg/Gauss]


function calc_env_chi(capsule_mass::Float64, tape_mass::Float64; T::Float64)::Float64
    A = 0.05275
    K1 = -2.0817e-8
    K2 = 1.25365e-8
    K3 = -1.07942e-12
    y = K1 + K2/T + K3*T
    return y * (capsule_mass + tape_mass) / A
end


function get_mag(dpath::String, mass_dict::Dict{String, Float64};
                corr_env::Bool=true,
                cols::Dict{String, String}=Dict(
                    "T"         =>"Temperature (K)",
                    "B"         =>"Field (Oe)",
                    "MAG"       =>"Long Moment (emu)",
                    "MAG_ERR"   =>"Long Scan Std Dev",
                    )
                )::DataFrame
    df_raw = parse_qdesign_file(dpath)

    molecular_weight    = mass_dict["molecular_weight"]
    sample_mass_g       = mass_dict["sample_mass_g"]
    mag_ions            = mass_dict["mag_ions"]
    pascal_corr         = mass_dict["pascal_corr"]
    capsule_mass_g      = mass_dict["capsule_mass_g"]
    tape_mass_g         = mass_dict["tape_mass_g"]

    temps       = df_raw[!, cols["T"]]
    field       = df_raw[!, cols["B"]]
    magmom      = df_raw[!, cols["MAG"]]    * (molecular_weight / (sample_mass_g * mag_ions))
    errmagmom   = df_raw[!, cols["MAG_ERR"]]* (molecular_weight / (sample_mass_g * mag_ions))

    if corr_env
        corr = field .* ([calc_env_chi(capsule_mass_g, tape_mass_g, T=t) for t in temps] .+ pascal_corr)
    else
        corr = 0.0
    end
    magmom .-= corr
    magmom .*= (1/(NA*muB))
    errmagmom .-= corr
    errmagmom .*= (1/(NA*muB))
    field .*= 1e-4
    df = DataFrame(T=temps, B=field, MAG=magmom, MAG_ERR=errmagmom)
    dropmissing!(df)
    sort!(df, :B)
    return df
end


function get_chi(dpath::String, mass_dict::Dict{String, Float64};
                corr_env::Bool=true,
                cols::Dict{String, String}=Dict(
                    "T"         =>"Temperature (K)",
                    "B"         =>"Field (Oe)",
                    "MAG"       =>"Long Moment (emu)",
                    "MAG_ERR"   =>"Long Scan Std Dev",
                    )
                )::DataFrame
    df_raw = parse_qdesign_file(dpath)

    molecular_weight    = mass_dict["molecular_weight"]
    sample_mass_g       = mass_dict["sample_mass_g"]
    mag_ions            = mass_dict["mag_ions"]
    pascal_corr         = mass_dict["pascal_corr"]
    capsule_mass_g      = mass_dict["capsule_mass_g"]
    tape_mass_g         = mass_dict["tape_mass_g"]

    temps   = df_raw[!, cols["T"]]
    field   = df_raw[!, cols["B"]]
    chi     = df_raw[!, cols["MAG"]]    * (molecular_weight / (sample_mass_g * mag_ions)) ./ field
    errchi  = df_raw[!, cols["MAG_ERR"]]* (molecular_weight / (sample_mass_g * mag_ions)) ./ field

    if corr_env
        corr = [calc_env_chi(capsule_mass_g, tape_mass_g, T=t) for t in temps] .+ pascal_corr
    else
        corr = 0.0
    end
    chi .-= corr
    errchi .-= corr
    field .*= 1e-4
    df = DataFrame(T=temps, B=field, CHI=chi, CHI_ERR=errchi)
    dropmissing!(df)
    sort!(df, :T)
    return df
end


function get_hc(dpath::String, mass_dict::Dict{String, Float64};
                cols::Dict{String, String}=Dict(
                    "T"         =>"Sample Temp (Kelvin)",
                    "B"         =>"Field (Oersted)",
                    "HC"        =>"Total HC (\xb5J/K)",
                    "HC_ERR"    =>"Total HC Err (\xb5J/K)",
                    )
                )::DataFrame
    df_raw = parse_qdesign_file(dpath)

    molecular_weight    = mass_dict["molecular_weight"] # g/mol
    sample_mass_g       = mass_dict["sample_mass_g"]    # g

    temps   = df_raw[!, cols["T"]]
    field   = df_raw[!, cols["B"]]
    hc      = df_raw[!, cols["HC"]]     * 1e-6 * (molecular_weight/sample_mass_g) # microJ/K ---> J/mol/K
    hc_err  = df_raw[!, cols["HC_ERR"]] * 1e-6 * (molecular_weight/sample_mass_g) # microJ/K ---> J/mol/K

    df = DataFrame(T=temps, B=field, HC=hc, HC_ERR=hc_err)
    dropmissing!(df)
    sort!(df, :T)
    return df
end


function make_massdict()::Dict{String, Float64}
    mdict = Dict(
        "molecular_weight"  => 0.0,
        "sample_mass_g"     => 0.0,
        "mag_ions"          => 0.0,
        "pascal_corr"       => 0.0,
        "capsule_mass_g"    => 0.0,
        "tape_mass_g"       => 0.0,
    )
    return mdict
end