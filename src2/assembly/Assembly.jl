module MolecularEvolution

using Base: max

struct FloydWarshall
    dist::Matrix{Int}
    next::Matrix{Int}
end

n(fw::FloydWarshall) = size(fw.dist, 1)

(fw::FloydWarshall)(i::Int, j::Int) = fw.dist[i, j]

"""
floydWarshallAlgorithm(x)

x is an adjacency list with fields:
- v::Vector{Int} : source vertices (0-based)
- w::Vector{Int} : destination vertices (0-based)
- label::Vector{Int} : edge weights

Returns a FloydWarshall struct with dist and next matrices.
    """
    function floydWarshallAlgorithm(x)
        nverts = max(maximum(x.v), maximum(x.w)) + 1
        max_weight = sum(x.label)

        dist = fill(max_weight, nverts, nverts)
        next = fill(-1, nverts, nverts)

        for i in eachindex(x.v)
            u = x.v[i] + 1
            v = x.w[i] + 1
            w = x.label[i]
            dist[u, v] = w
            next[u, v] = v
        end

        for k in 1:nverts, i in 1:nverts, j in 1:nverts
            if dist[i, k] + dist[k, j] < dist[i, j]
                dist[i, j] = dist[i, k] + dist[k, j]
                next[i, j] = next[i, k]
            end
        end

        FloydWarshall(dist, next)
    end

    """
    floydWarshallPathReconstruction(u, v, fw)

    Reconstructs the path from u to v using FloydWarshall result fw.
    Returns a vector of node indices (1-based).
    """
    function floydWarshallPathReconstruction(u::Int, v::Int, fw::FloydWarshall)
        @assert 1 ≤ u ≤ n(fw) "u=$u out of bounds"
        @assert 1 ≤ v ≤ n(fw) "v=$v out of bounds"

        path = Int[]
        fw.next[u, v] == -1 && return path

        push!(path, u)
        current = u
        while current != v
            current = fw.next[current, v]
            push!(path, current)
        end

        path
    end

end # module
