#!/usr/bin/env julia

# Using AbstractString allows the struct to hold both String and SubString
struct SuffixTree
    text::AbstractString
    adj::Dict{Int, Vector{Tuple{Int, Int, Int}}}
end

"""
build_suffix_tree(text::AbstractString)
Constructs a Suffix Tree by first building a Trie of all suffixes
and then compressing paths.
"""
function build_suffix_tree(text::AbstractString)
    # 1. Build a Suffix Trie
    # adj[parent_node] = [(child_node, start_index, length), ...]
    adj = Dict{Int, Vector{Tuple{Int, Int, Int}}}()
    adj[1] = []
    next_node = 2

    n = length(text)
    for i in 1:n
        current_node = 1
        # Iterate through the suffix starting at i
        for j in i:n
            char = text[j]
            found = false
            for (child, start, len) in adj[current_node]
                if text[start] == char
                    current_node = child
                    found = true
                    break
                end
            end

            if !found
                adj[next_node] = []
                push!(adj[current_node], (next_node, j, 1))
                current_node = next_node
                next_node += 1
            end
        end
    end

    # 2. Compress the Trie into a Suffix Tree
    compressed_adj = Dict{Int, Vector{Tuple{Int, Int, Int}}}()

    function compress!(u)
        compressed_adj[u] = []
        for (v, start, len) in adj[u]
            curr_v, curr_start, curr_len = v, start, len
            # If a node has exactly one child, merge the edges
            while length(adj[curr_v]) == 1
                next_v, next_start, next_len = adj[curr_v][1]
                curr_v = next_v
                curr_len += next_len
            end
            push!(compressed_adj[u], (curr_v, curr_start, curr_len))
            compress!(curr_v)
        end
    end

    compress!(1)
    return SuffixTree(text, compressed_adj)
end

"""
get_edge_substrings(st::SuffixTree)
Extracts the actual substring from the original text for each edge.
    """
    function get_edge_substrings(st::SuffixTree)
        edge_strings = String[]
        for u in keys(st.adj)
            for (v, start, len) in st.adj[u]
                # Slicing the AbstractString and converting to String for output
                push!(edge_strings, String(st.text[start:start+len-1]))
            end
        end
        return edge_strings
    end

    function main()
        input_path = "data/stepic_9d.txt"
        if !isfile(input_path)
            println("File not found: $input_path")
            return
        end

        # strip() returns a SubString{String}
        text = strip(read(input_path, String))

        st = build_suffix_tree(text)
        edges = get_edge_substrings(st)

        # Output formatting
        output_text = join(edges, "\n")
        println(output_text)

        mkpath("output")
        write("output/Assignment_09D.txt", output_text)
    end

    # Execution check
    if abspath(PROGRAM_FILE) == @__FILE__
        main()
    end

