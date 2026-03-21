using DataStructures

# --- Helper: Hierholzer's Algorithm ---
function eulerian_cycle(adj)
    # Start at any node that has outgoing edges
    start_node = first(keys(adj))
    stack = [start_node]
    cycle = []
    
    # Work on a copy to preserve original dictionary
    local_adj = Dict(k => copy(v) for (k, v) in adj)

    while !isempty(stack)
        u = stack[end]
        if haskey(local_adj, u) && !isempty(local_adj[u])
            v = popfirst!(local_adj[u]) 
            push!(stack, v)
        else
            push!(cycle, pop!(stack))
        end
    end
    return reverse(cycle)
end

# --- Eulerian Path Logic ---
function find_eulerian_path(adj)
    in_degree = DefaultDict{Any, Int}(0)
    out_degree = DefaultDict{Any, Int}(0)
    nodes = Set()

    for (u, targets) in adj
        push!(nodes, u)
        out_degree[u] += length(targets)
        for v in targets
            push!(nodes, v)
            in_degree[v] += 1
        end
    end

    start_node, end_node = nothing, nothing
    for node in nodes
        if out_degree[node] > in_degree[node]
            start_node = node
        elseif in_degree[node] > out_degree[node]
            end_node = node
        end
    end

    if !isnothing(start_node) && !isnothing(end_node)
        if !haskey(adj, end_node) adj[end_node] = [] end
        push!(adj[end_node], start_node)
    end

    cycle = eulerian_cycle(adj)

    if !isnothing(start_node)
        for i in 1:(length(cycle)-1)
            if cycle[i] == end_node && cycle[i+1] == start_node
                return [cycle[i+1:end-1]; cycle[1:i]]
            end
        end
    end
    return cycle
end

# --- Main Solve Function ---
function solve()
    input_path = "data/stepic_5d.txt"
    if !isfile(input_path)
        println("Error: File not found at $input_path")
        return
    end

    lines = readlines(input_path)
    if isempty(lines)
        println("Error: File is empty.")
        return
    end

    # 1. Improved Parser: Find k and d anywhere in the first few lines
    metadata = Int[]
    paired_reads = Vector{Vector{String}}()

    for line in lines
        line = strip(line)
        if isempty(line) continue end
        
        if contains(line, "|")
            push!(paired_reads, split(line, "|"))
        else
            # Extract all integers found in non-DNA lines
            parts = split(line)
            for p in parts
                val = tryparse(Int, p)
                if val !== nothing
                    push!(metadata, val)
                end
            end
        end
    end

    if length(metadata) < 2
        # FALLBACK: If k and d aren't in the file, we infer k from the first read
        if !isempty(paired_reads)
            k = length(paired_reads[1][1])
            d = 100 # Default d if missing; change this if you know your d
            @warn "k and d not found in header. Inferred k=$k, using d=$d."
        else
            error("Could not find metadata OR DNA reads in the file.")
        end
    else
        k, d = metadata[1], metadata[2]
    end

    # 2. Build De Bruijn Graph
    adj = Dict{Tuple{String, String}, Vector{Tuple{String, String}}}()
    for pair in paired_reads
        # Use explicit k value for prefix/suffix
        u = (pair[1][1:k-1], pair[2][1:k-1])
        v = (pair[1][2:k],   pair[2][2:k])
        
        if !haskey(adj, u) adj[u] = [] end
        push!(adj[u], v)
    end

    # 3. Find Path
    path = find_eulerian_path(adj)

    # 4. String Reconstruction
    string1 = path[1][1] * join([node[1][end] for node in path[2:end]])
    string2 = path[1][2] * join([node[2][end] for node in path[2:end]])

    # The overlap logic: the prefix of string2 aligns with string1 after k+d positions
    result = string1 * string2[end-k-d+1:end]

    println("Success! Result length: ", length(result))
    println(result)
    
    if !isdir("output") mkpath("output") end
    write("output/Assignment_05D.txt", result)
end

solve()

