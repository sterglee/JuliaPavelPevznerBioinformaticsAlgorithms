using DataStructures

# -----------------------------
# Helper: Eulerian Cycle
# -----------------------------
function eulerian_cycle(graph::Dict{Int, Vector{Int}})
    # Make a deep copy so we can mutate safely
    graph = deepcopy(graph)

    # Stack for current path and final cycle
    stack = []
    cycle = []

    # Start from any node with edges
    start_node = first(keys(graph))
    push!(stack, start_node)

    while !isempty(stack)
        v = last(stack)

        if haskey(graph, v) && !isempty(graph[v])
            u = pop!(graph[v])   # remove edge v -> u
            push!(stack, u)
        else
            push!(cycle, pop!(stack))
        end
    end

    return reverse(cycle)
end

# -----------------------------
# Main Execution
# -----------------------------
    input_file = "data/stepic_4e.txt"
    output_dir = "output"
    output_file = joinpath(output_dir, "Assignment_04E.txt")

    # Ensure directories exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Robust parsing
    edges = Dict{Int, Vector{Int}}()

    if isfile(input_file)
        open(input_file) do input_data
            for line in eachline(input_data)
                line = strip(line)
                if isempty(line)
                    continue
                end

                # Split the "Node -> Targets" format
                parts = split(line, " -> ")

                u = parse(Int, parts[1])

                # Parse comma-separated targets into Vector{Int}
                v_list = parse.(Int, split(parts[2], ","))

                edges[u] = v_list
            end
        end

        # Find the cycle
        path = eulerian_cycle(edges)

        # Join with "->" for final output
        result_str = join(path, "->")

        println("Success: Cycle found.")
        println(result_str)

        open(output_file, "w") do f
            write(f, result_str)
        end
    else
        println("Error: Input file $input_file not found.")
    end

