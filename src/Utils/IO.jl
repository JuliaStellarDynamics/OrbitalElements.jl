using HDF5 # To have access to .hf5 formats


"""
write all the parameters to a file
"""
function WriteParameters(filename::String,
                         Parameters::OrbitsParameters,
                         mode::String="r+")

    h5open(filename, mode) do file
        WriteParameters(file,Parameters)
    end
end

function WriteParameters(file::HDF5.File,
                         Parameters::OrbitsParameters)

    group = create_group(file,"OrbitsParameters")
    for i = 1:fieldcount(OrbitsParameters)
        varname = string(fieldname(OrbitsParameters,i))
        try write(group,varname,getfield(Parameters,i)) catch; println("Unable to write parameter: "*varname) end
    end
end