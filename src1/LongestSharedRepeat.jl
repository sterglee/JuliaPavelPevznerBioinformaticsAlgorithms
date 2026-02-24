#!/usr/bin/env julia

"""
Node
Represents a node in the Generalized Suffix Tree.
- `children`: maps a starting character to an edge.
- `word_indices`: keeps track of which input strings (1 or 2) pass through this node.
"""
mutable struct Node
    children::Dict{Char, Tuple{Int, Int, Int}} # Char => (child_node_id, start_idx, len)
    word_sources::Set{Int}
    depth::Int
end

Node(d=0) = Node(Dict(), Set{Int}(), d)

"""
build_generalized_suffix_tree(strings::Vector{<:AbstractString})
Builds a GST by inserting suffixes from all provided strings.
To distinguish strings, we conceptually treat them as concatenated with unique separators.
"""
function build_generalized_suffix_tree(strings)
    root = Node()
    nodes = [root]

    for (s_idx, original_s) in enumerate(strings)
        # We append a unique terminal for each string to handle the "Generalized" aspect
        s = original_s * (s_idx == 1 ? "#" : "\$")
        n = length(s)

        for i in 1:n
            current_node_idx = 1
            suffix = @view s[i:end]

            for char in suffix
                # Mark that this node is visited by string s_idx
                push!(nodes[current_node_idx].word_sources, s_idx)

                if haskey(nodes[current_node_idx].children, char)
                    child_idx, start, len = nodes[current_node_idx].children[char]
                    current_node_idx = child_idx
                else
                    # Create new node
                    new_node = Node(nodes[current_node_idx].depth + 1)
                    push!(nodes, new_node)
                    new_idx = length(nodes)

                    # For a simple Trie-based build, len is 1.
                    # We will compress later or just use the depth for LCS.
                    nodes[current_node_idx].children[char] = (new_idx, i, 1)
                    current_node_idx = new_idx
                end
            end
            # Mark the leaf
            push!(nodes[current_node_idx].word_sources, s_idx)
        end
    end
    return nodes, strings
end

"""
find_lcs(nodes, strings)
Finds the longest string path that is shared by all input word sources.
"""
function find_lcs(nodes, strings)
    num_strings = length(strings)
    best_depth = -1
    best_path = ""

    # Helper to reconstruct path from root to node (simplified for this Trie approach)
    # In a compressed tree, we'd track parent pointers.
    # Here, we'll perform a DFS to find the deepest node with all word_sources.

    best_node_path = ""

    function dfs(u_idx, current_path)
        u = nodes[u_idx]

        # If this node isn't shared by all strings, stop
        if length(u.word_sources) < num_strings
            return
        end

        # Check if this is the deepest shared node found so far
        # We exclude the terminal symbols from the length
        clean_path = replace(current_path, r"#|\$" => "")
        if length(clean_path) > best_depth
            best_depth = length(clean_path)
            best_path = clean_path
        end

        for (char, (v_idx, start, len)) in u.children
            dfs(v_idx, current_path * char)
        end
    end

    dfs(1, "")
    return best_path
end

function main()
    input_path = "data/stepic_9e.txt"
    if !isfile(input_path)
        return
    end

    lines = filter(!isempty, strip.(readlines(input_path)))
    if length(lines) < 2
        return
    end

    # 1. Build the tree
    nodes, original_strings = build_generalized_suffix_tree(lines)

    # 2. Find Longest Common Substring
    result = find_lcs(nodes, original_strings)

    # 3. Output
    println(result)
    mkpath("output")
    write("output/Assignment_09E.txt", result)
end

    main()

