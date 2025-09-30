
"""
    WriteParameters(filename,params,mode)
    WriteParameters(file,params)

write all the parameters to a file
"""
function WriteParameters(filename::String,
                         params::OrbitalParameters=OrbitalParameters(),
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,params)
    end
end

function WriteParameters(file::HDF5.File,
                         params::OrbitalParameters=OrbitalParameters())

    group = create_group(file,"OrbitalParameters")
    for i = 1:fieldcount(OrbitalParameters)
        varname = string(fieldname(OrbitalParameters,i))
        if !isascii(varname)
            varname = NonAsciiHandle(varname)
        end
        try write(group,varname,getfield(params,i)) catch; println("Unable to write parameter: "*varname) end
    end
end

"""
    NonAsciiHandle(x)

convert some extra unicode characters to ascii
"""
function NonAsciiHandle(x::String)::String

    return replace(x,"Ω"=>"Omega","₀"=>"0","α"=>"alpha","ε"=>"eps",!isascii=>"")
end