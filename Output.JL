module Output

function writeStringToFile(s::String, path::String)
    open(path, "w") do io
        write(io, s)
    end
    nothing
end

end # module
