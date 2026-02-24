module Assembly

using BioSequences

"""
composition(k::Int, text::AbstractString) -> Vector{String}

Solve the String Composition Problem.

Input: An integer k and a string Text.
Output: Compositionk(Text) (the k-mers can be provided in any order).

Sample Input:
5
CAATCCAAC
Sample Output:
CAATC
AATCC
ATCCA
TCCAA
CCAAC
"""
function composition(k::Int, text::AbstractString)
    kmers = [text[i:i+k-1] for i in 1:(length(text)-k+1)]
        unique_sorted_kmers = sort(unique(kmers))
        return unique_sorted_kmers
    end

    """
    genomePath(pathElems::Vector{String}) -> String

    Solve the String Spelled by a Genome Path Problem.

    genomePath(["ACCGA","CCGAA", "CGAGC", "AAGCT"])

    Sample Input:
    ACCGA
    CCGAA
    CGAAG
    GAAGC
    AAGCT
    Sample Output:
    ACCGAAGCT
    """
    function genomePath(pathElems::Vector{String})
        @assert length(pathElems) > 1
        first = pathElems[1]
        @assert length(first) > 0

        seq = Vector{Char}(undef, length(first) + length(pathElems) - 1)
        for i in 1:length(seq)
            if i <= length(first)
                seq[i] = first[i]
            else
                seq[i] = last(pathElems[i - length(first) + 1])
            end
        end
        return String(seq)
    end

    """
    overlap(patterns::Vector{String}) -> Vector{Tuple{String,String}}
overlap(["ATGCG","GCATG", "CATGC", "AGGCA","GGCAT"])
    Solve the Overlap Graph Problem.

    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns), in the form of an adjacency list.

    Sample Input:
    ATGCG
    GCATG
    CATGC
    AGGCA
    GGCAT
    Sample Output:
    CATGC -> ATGCG
    GCATG -> CATGC
    GGCAT -> GCATG
    AGGCA -> GGCAT
    """
    function prefixes(patterns::Vector{String})
        prefix_map = Dict{String, Vector{Int}}()
        for (i, pattern) in enumerate(patterns)
            prefix = pattern[1:end-1]
            push!(get!(prefix_map, prefix, Int[]), i)
        end
        return prefix_map
    end

    function overlap(patterns::Vector{String})
        prefix_map = prefixes(patterns)
        unique_sorted_patterns = sort(unique(patterns))
        edges = Vector{Tuple{String,String}}()
        for pattern in unique_sorted_patterns
            suffix = pattern[2:end]
            if haskey(prefix_map, suffix)
                for i in prefix_map[suffix]
                    push!(edges, (pattern, patterns[i]))
                end
            end
        end
        return edges
    end

    """
    deBrujinGraph(k::Int, text::AbstractString) -> Vector{Tuple{String, Vector{String}}}

    Solve the De Bruijn Graph from a String Problem.

    Input: An integer k and a string Text.
    Output: DeBruijnk(Text), in the form of an adjacency list.

    Sample Input:
    4
    AAGATTCTCTAAGA
    Sample Output:
    AAG -> AGA,AGA
    AGA -> GAT
    ATT -> TTC
    CTA -> TAA
    CTC -> TCT
    GAT -> ATT
    TAA -> AAG
    TCT -> CTA,CTC
    TTC -> TCT
    """
    function getKmers(text::AbstractString, k::Int)
        return [text[i:i+k-1] for i in 1:(length(text)-k+1)]
        end

        function deBrujinGraphFromKmers(kmers::Vector{String})
            adjacency = Dict{String, Vector{String}}()
            for kmer in kmers
                prefix = kmer[1:end-1]
                suffix = kmer[2:end]
                push!(get!(adjacency, prefix, String[]), suffix)
            end
            # Convert to sorted vector of tuples
            result = [(k, adjacency[k]) for k in sort(collect(keys(adjacency)))]
                return result
            end

            function deBrujinGraph(k::Int, text::AbstractString)
                kmers = getKmers(text, k)
                return deBrujinGraphFromKmers(kmers)
            end

        end # module

