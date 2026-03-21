using Printf
using DelimitedFiles
using DataStructures

# -----------------------------
# Trie Implementation
# -----------------------------
mutable struct TrieNode
    children::Dict{Char, TrieNode}
    id::Int
    is_root::Bool

    function TrieNode(id::Int, is_root::Bool=false)
        new(Dict{Char, TrieNode}(), id, is_root)
    end
end

mutable struct Trie
    root::TrieNode
    edges::Vector{Tuple{Tuple{Int, Int}, Char}}
    node_count::Int

    function Trie()
        root = TrieNode(1, true)
        new(root, [], 1)
    end
end

function insert!(t::Trie, word::AbstractString)
    current = t.root

    for c in word
        if !haskey(current.children, c)
            t.node_count += 1
            new_node = TrieNode(t.node_count)
            current.children[c] = new_node

            # Store edge
            push!(t.edges, ((current.id, new_node.id), c))
        end
        current = current.children[c]
    end
end

function Trie(words::AbstractVector{<:AbstractString})
    t = Trie()
    for w in words
        insert!(t, w)
    end
    return t
end

# -----------------------------
# Core Function (FIXED)
# -----------------------------
function trie_edges(words::AbstractVector{<:AbstractString})
    """Returns the edges of a trie constructed from the given words."""

    t = Trie(words)

    adjacency_format(edge) = @sprintf("%d %d %c", edge[1][1], edge[1][2], edge[2])

    return map(adjacency_format, t.edges)
end

# -----------------------------
# Main Execution
# -----------------------------
function main()
    input_file = "data/stepic_9a.txt"
    output_dir = "output"
    output_file = joinpath(output_dir, "Assignment_09A.txt")

    if !isdir(output_dir)
        mkpath(output_dir)
    end

    words = readlines(input_file) .|> strip

    adjacency_list = trie_edges(words)

    println(join(adjacency_list, "\n"))

    open(output_file, "w") do io
        write(io, join(adjacency_list, "\n"))
    end
end

main()

