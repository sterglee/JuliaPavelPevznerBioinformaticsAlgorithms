#!/usr/bin/env julia

"""
Node
- `children`: Char => (child_node_id, start_index, length)
- `sources`: Set of string indices (1 for String A, 2 for String B) that pass through this node.
    - `parent`: ID of the parent node.
    - `edge_info`: Tuple of (start_idx, len) to reconstruct the edge label.
    """
    mutable struct Node
        children::Dict{Char, Tuple{Int, Int, Int}}
        sources::Set{Int}
        parent::Int
        edge_info::Tuple{Int, Int}
        depth::Int
    end

    Node(p=0, info=(0,0), d=0) = Node(Dict(), Set{Int}(), p, info, d)

    """
    build_gst(s1::AbstractString, s2::AbstractString)
    Constructs a simple Generalized Suffix Tree for two strings.
        """
        function build_gst(s1, s2)
            # Append unique terminals to ensure suffixes are explicit
            # String 1 ends in '#', String 2 ends in '$'
            strings = [s1 * "#", s2 * "\$"]
            nodes = [Node()] # Root at index 1

            for (idx, s) in enumerate(strings)
                n = length(s)
                for i in 1:n
                    curr = 1
                    suffix = @view s[i:end]
                    for (char_idx, char) in enumerate(suffix)
                        push!(nodes[curr].sources, idx)

                        if haskey(nodes[curr].children, char)
                            curr, _, _ = nodes[curr].children[char]
                        else
                            # Create new node
                            new_node = Node(curr, (i + char_idx - 1, 1), nodes[curr].depth + 1)
                            push!(nodes, new_node)
                            new_idx = length(nodes)
                            nodes[curr].children[char] = (new_idx, i + char_idx - 1, 1)
                            curr = new_idx
                        end
                    end
                    push!(nodes[curr].sources, idx)
                end
            end
            return nodes, strings
        end

        """
        solve_shortest_unique(nodes, strings)
        Finds the shortest substring in strings[1] that does not appear in strings[2].
        """
        function solve_shortest_unique(nodes, strings)
            s1_full = strings[1] * "#"
            shortest_len = Inf
            result = ""

            # We look for nodes that are owned ONLY by source 1
            for (i, node) in enumerate(nodes)
                if node.sources == Set([1])
                    # The parent of this node MUST be shared or also unique.
                    # A unique substring is formed by taking the path to the parent
                    # (which is shared) plus the first character of the edge to this unique node.
                    parent_node = nodes[node.parent]

                    # Reconstruct path to parent
                    path_to_parent = ""
                    curr_p = node.parent
                    while curr_p > 1
                        p_node = nodes[curr_p]
                        start, len = p_node.edge_info
                        # We use the original string 1 to reconstruct
                        path_to_parent = s1_full[start:start+len-1] * path_to_parent
                        curr_p = p_node.parent
                    end

                    # Add the first character of the edge that made it unique
                    edge_start = node.edge_info[1]
                    first_char_of_edge = s1_full[edge_start]

                    # Skip the terminal '#' character if it's the unique part
                    if first_char_of_edge == '#' continue end

                    candidate = path_to_parent * first_char_of_edge
                    if length(candidate) < shortest_len
                        shortest_len = length(candidate)
                        result = candidate
                    end
                end
            end
            return result
        end

        function main()
            input_path = "data/stepic_9f.txt"
            if !isfile(input_path) return end

            lines = filter(!isempty, strip.(readlines(input_path)))
            if length(lines) < 2 return end

            # Build Tree
            nodes, strings = build_gst(lines[1], lines[2])

            # Logic: Shortest string in A not in B
            ans = solve_shortest_unique(nodes, lines)

            println(ans)
            mkpath("output")
            write("output/Assignment_09F.txt", ans)
        end

        if abspath(PROGRAM_FILE) == @__FILE__
            main()
        end

