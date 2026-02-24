using DelimitedFiles

"""
Solves the Longest Path in a DAG problem.
"""
function solve_longest_path()
    # 1. Setup File Paths
    input_path = "data/stepic_6d.txt"
    output_path = "output/Assignment_06D.txt"

    if !isfile(input_path)
        println("Error: Input file not found at $input_path")
        return
    end

    # 2. Parse Input
    lines = readlines(input_path)
    source = parse(Int, lines[1])
    sink = parse(Int, lines[2])

    adj = Dict{Int, Vector{Tuple{Int, Int}}}()
    in_degree = Dict{Int, Int}()
    all_nodes = Set{Int}([source, sink])

    for i in 3:length(lines)
        line = strip(lines[i])
        isempty(line) && continue

        parts = split(line, "->")
        u = parse(Int, parts[1])
        v_parts = split(parts[2], ":")
        v = parse(Int, v_parts[1])
        w = parse(Int, v_parts[2])

        push!(all_nodes, u, v)
        push!(get!(adj, u, Tuple{Int, Int}[]), (v, w))
        in_degree[v] = get(in_degree, v, 0) + 1
        get!(in_degree, u, 0)
    end

    # 3. Topological Sort (Kahn's Algorithm)

    queue = [n for n in all_nodes if get(in_degree, n, 0) == 0]
        top_order = Int[]
        temp_in_degree = copy(in_degree)

        while !isempty(queue)
            u = popfirst!(queue)
            push!(top_order, u)
            if haskey(adj, u)
                for (v, _) in adj[u]
                    temp_in_degree[v] -= 1
                    if temp_in_degree[v] == 0
                        push!(queue, v)
                    end
                end
            end
        end

        # 4. Longest Path Calculation (DP)
        # Using a very large negative number for unreachable nodes
        dist = Dict{Int, Int64}(n => -10^15 for n in all_nodes)
            backtrack = Dict{Int, Int}()
            dist[source] = 0

            for u in top_order
                if dist[u] <= -10^14
                    continue
                end

                if haskey(adj, u)
                    for (v, w) in adj[u]
                        if dist[u] + w > dist[v]
                            dist[v] = dist[u] + w
                            backtrack[v] = u
                        end
                    end
                end
            end

            # 5. Path Reconstruction and Output
            # Check if sink is reachable
            if !haskey(dist, sink) || dist[sink] <= -10^14
                println("No path exists from $source to $sink.")
                return
            end

            path = [sink]
            while path[1] != source
                pushfirst!(path, backtrack[path[1]])
            end

            # Define variables inside the function to avoid "not defined in Main" errors
            max_dist = dist[sink]
            path_str = join(path, "->")
            final_result = "$max_dist\n$path_str"

            # Print and Write
            println(final_result)

            mkpath(dirname(output_path))
            write(output_path, final_result)
            println("\nSuccess: Result saved to $output_path")
        end

        # Execute the function
        solve_longest_path()

