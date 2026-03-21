using DataStructures

"""
Finds an Eulerian Cycle in a directed graph using Hierholzer's Algorithm.
"""
function eulerian_cycle(graph::Dict{String, Vector{String}})
    # Start at any node that has outgoing edges
    start_node = first(keys(graph))
    stack = [start_node]
    circuit = String[]

    # We make a local copy of the adjacency list to avoid destroying the input
    adj = Dict(k => copy(v) for (k, v) in graph)

    while !isempty(stack)
        u = stack[end]
        if haskey(adj, u) && !isempty(adj[u])
            # Move to the next node and keep the edge on the stack
            v = popfirst!(adj[u]) 
            push!(stack, v)
        else
            # Backtrack: add to circuit when no more outgoing edges exist
            push!(circuit, pop!(stack))
        end
    end
    
    # Hierholzer's returns the path in reverse
    return reverse(circuit)
end

"""
Generates a binary string of length 2^k that contains all 2^k binary k-mers.
"""
function universal_string_problem(k)
    if k == 1 return "01" end
    
    # 1. Construct the De Bruijn Graph
    # Nodes are (k-1)-mers, Edges are k-mers
    universal_dict = Dict{String, Vector{String}}()
    
    # Generate all possible binary strings of length k
    # Iterators.product('0':'1', ...) generates tuples of characters
    for bits in Iterators.product(fill('0':'1', k)...)
        kmer = join(bits)
        prefix = kmer[1:end-1]
        suffix = kmer[2:end]
        
        if !haskey(universal_dict, prefix)
            universal_dict[prefix] = String[]
        end
        push!(universal_dict[prefix], suffix)
    end
    
    # 2. Find Eulerian Cycle
    cycle = eulerian_cycle(universal_dict)
    
    # 3. Form the string
    # In a circular string, each node in the cycle provides one character.
    # We take the first character of the first 2^k nodes.
    return join(node[1] for node in cycle[1:end-1])
end

# Example usage:
k = 4
result = universal_string_problem(k)
println("k=$k Universal String: ", result)
println("Length: ", length(result)) # Should be 2^k = 16

