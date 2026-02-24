module MolecularEvolution

using Base: max

struct FloydWarshall
    dist::Vector{Vector{Int}}
    next::Vector{Vector{Int}}
    function FloydWarshall(dist::Vector{Vector{Int}}, next::Vector{Vector{Int}})
        new(dist, next)
    end
end

function n(fw::FloydWarshall)
    return length(fw.dist)
end

function (fw::FloydWarshall)(i::Int, j::Int)
    return fw.dist[i][j]
end

"""
floydWarshallAlgorithm(x)

x is an adjacency list with fields:
- v::Vector{Int} : source vertices
- w::Vector{Int} : destination vertices
- label::Vector{Int} : weights of edges

Returns a FloydWarshall struct with dist and next matrices.
    """
    function floydWarshallAlgorithm(x)
        n = max(maximum(x.v), maximum(x.w)) + 1
        max_weight = sum(x.label)
        dist = [max_weight for _ in 1:n, _ in 1:n]
            next = [-1 for _ in 1:n, _ in 1:n]

                for i in eachindex(x.v)
                    u = x.v[i] + 1  # Julia is 1-based indexing
                    v = x.w[i] + 1
                    w = x.label[i]
                    dist[u][v] = w
                    next[u][v] = v
                end

                for k in 1:n
                    for i in 1:n
                        for j in 1:n
                            if dist[i][k] + dist[k][j] < dist[i][j]
                                dist[i][j] = dist[i][k] + dist[k][j]
                                next[i][j] = next[i][k]
                            end
                        end
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
                @assert 1 <= u <= n(fw) "u=$u out of bounds"
                @assert 1 <= v <= n(fw) "v=$v out of bounds"
                path = Int[]
                if fw.next[u][v] >= 0
                    push!(path, u)
                    currentNode = u
                    while currentNode != v
                        currentNode = fw.next[currentNode][v]
                        push!(path, currentNode)
                    end
                end
                return path
            end

        end # module
