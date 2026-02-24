#!/usr/bin/env julia

# Fast FASTA parser (single or multi-record)
function read_fasta_sequence(path::String)::String
    seq_buffer = IOBuffer()

    open(path, "r") do io
        for line in eachline(io)
            isempty(line) && continue
            if line[1] == '>'
                continue   # skip FASTA header
            end
            write(seq_buffer, strip(line))
        end
    end

    return String(take!(seq_buffer))
end


function most_frequent_kmers(dna::String, k::Int)
    n = length(dna)
    counts = Dict{SubString{String}, Int}()

    @inbounds for i in 1:(n - k + 1)
        kmer = SubString(dna, i, i + k - 1)
        counts[kmer] = get(counts, kmer, 0) + 1
    end

    max_count = maximum(values(counts))

    # Convert back to String only for output
    return [String(kmer) for (kmer, c) in counts if c == max_count]
    end


    # ---------------------------
    # Main
    # ---------------------------

    input_fasta = "data/sequence.fasta"
    k = 9  # change as needed

    dna = read_fasta_sequence(input_fasta)

    kmers = most_frequent_kmers(dna, k)

    result = join(kmers, " ")
    println(result)

    open("output/createDistanceMatrix.res", "w") do io
        write(io, result)
    end
