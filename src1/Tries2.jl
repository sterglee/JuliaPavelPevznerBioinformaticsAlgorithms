using Printf
using DataStructures

# -----------------------------
# Trie Implementation
# -----------------------------
mutable struct TrieNode
    children::Dict{Char, TrieNode}
    is_end::Bool

    function TrieNode()
        new(Dict{Char, TrieNode}(), false)
    end
end

mutable struct Trie
    root::TrieNode
end

function Trie(patterns::AbstractVector{<:AbstractString})
    t = Trie(TrieNode())
    for p in patterns
        insert!(t, p)
    end
    return t
end

function insert!(t::Trie, word::AbstractString)
    current = t.root
    for c in word
        if !haskey(current.children, c)
            current.children[c] = TrieNode()
        end
        current = current.children[c]
    end
    current.is_end = true
end

# -----------------------------
# Prefix Matching
# -----------------------------
function prefix_in_trie(t::Trie, text::AbstractString)
    current = t.root
    for c in text
        if current.is_end
            return true
        end
        if !haskey(current.children, c)
            return false
        end
        current = current.children[c]
    end
    return current.is_end
end

# -----------------------------
# Pattern Matching
# -----------------------------
function trie_pattern_matching(word::AbstractString, patterns::AbstractVector{<:AbstractString})
    """Returns the starting index of all locations in word where a string in patterns is a substring."""

    t = Trie(patterns)

    min_len = minimum(length.(patterns))
    n = length(word)

    indices = Int[]

    for i in 1:(n - min_len + 1)
        if prefix_in_trie(t, word[i:end])
            push!(indices, i - 1)  # 0-based indexing (Bioinformatics standard)
        end
    end

    return indices
end

# -----------------------------
# Main Execution
# -----------------------------
function main()
    """Main call. Reads, runs, and saves problem specific data."""

    input_file = "data/stepic_9b.txt"
    output_dir = "output"
    output_file = joinpath(output_dir, "Assignment_09B.txt")

    if !isdir(output_dir)
        mkpath(output_dir)
    end

    input_data = readlines(input_file)

    word = strip(input_data[1])
    patterns = strip.(input_data[2:end])

    pattern_indices = trie_pattern_matching(word, patterns)

    println(join(pattern_indices, " "))

    open(output_file, "w") do io
        write(io, join(pattern_indices, " "))
    end
end

# -----------------------------
# Suffix Tree (Simplified)
# -----------------------------
function suffix_tree_edges(word::AbstractString)
    """Returns edge substrings of a naive suffix tree (simplified)."""

    suffixes = [word[i:end] for i in 1:length(word)]

    edges = Set{String}()

    for s in suffixes
        for i in 1:length(s)
            push!(edges, s[1:i])
        end
    end

    return collect(edges)
end

# Run
main()

