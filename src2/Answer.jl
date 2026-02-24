module Answer

# Assuming these types/functions exist elsewhere, mirroring the Scala code:
# - AdjacencyList with fields: v, w, label
# - SuffixTree with fields: text, edges (edges has v, pos, len)
# - int2label(::Int)

export printVals, printSeq, printTuples, printDeBrujinGraph,
       printMatrix, printAdjacencyList, printTreeAdjacencyList,
       printSuffixTreeEdgeLabels

function printVals(x::AbstractVector{Int})
    io = IOBuffer()
    for s in x
        print(io, s, " ")
    end
    String(take!(io))
end

function printSeq(x::AbstractVector{String})
    io = IOBuffer()
    for s in x
        print(io, s, " ")
    end
    String(take!(io))
end

function printTuples(x::AbstractVector{Tuple{String,String}})
    io = IOBuffer()
    for i in 1:length(x)-1
        a, b = x[i]
        println(io, a, " -> ", b)
    end
    a, b = x[end]
    print(io, a, " -> ", b)
    String(take!(io))
end

function printDeBrujinGraph(x::AbstractVector{Tuple{String,AbstractVector{String}}})
    io = IOBuffer()
    n = length(x)
    for (i, (node, patterns)) in enumerate(x)
        print(io, node, " -> ")
        for (j, p) in enumerate(patterns)
            j > 1 && print(io, ",")
            print(io, p)
        end
        i < n && print(io, "\n")
    end
    String(take!(io))
end

function printMatrix(x::AbstractVector{<:AbstractVector{Int}}, sep::String = "\t")
    io = IOBuffer()
    for i in eachindex(x)
        row = x[i]
        for j in eachindex(row)
            print(io, row[j])
            j != lastindex(row) && print(io, sep)
        end
        i != lastindex(x) && print(io, "\n")
    end
    String(take!(io))
end

function printAdjacencyList(x)
    io = IOBuffer()
    indices = sort(eachindex(x.v), by = i -> x.v[i])
    for i in indices
        println(io,
            x.v[i], "->", x.w[i], ":",
            int2label(x.label[i])
        )
    end
    str = String(take!(io))
    chop(str; tail = 1)  # remove final newline
end

function printTreeAdjacencyList(x)
    io = IOBuffer()
    indices = sort(eachindex(x.v), by = i -> x.v[i])
    for i in indices
        println(io,
            x.v[i], "->", x.w[i], ":",
            x.label[i]
        )
    end
    str = String(take!(io))
    chop(str; tail = 1)
end

function printSuffixTreeEdgeLabels(x)
    io = IOBuffer()
    for i in eachindex(x.edges.v)
        start = x.edges.pos[i] + 1
        stop  = start + x.edges.len[i] - 1
        println(io, x.text[start:stop])
    end
    str = String(take!(io))
    chop(str; tail = 1)
end

end # module
