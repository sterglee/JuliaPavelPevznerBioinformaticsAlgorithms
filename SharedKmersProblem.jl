#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Problem Title: Shared k-mers Problem
=#

"""
reverse_complement(s::AbstractString)

Returns the reverse complement of a DNA string.
"""
function reverse_complement(s::AbstractString)
    complement = Dict('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C',
                      'a'=>'t', 't'=>'a', 'c'=>'g', 'g'=>'c')
    # reverse(s) and then map the complement
    return map(c -> complement[c], reverse(s))
end

"""
shared_kmers(k::Int, dna1::String, dna2::String)

Returns a set of (i, j) tuples where dna1[i:i+k-1] matches dna2[j:j+k-1]
or its reverse complement.
"""
function shared_kmers(k, dna1, dna2)
    # Using a Dict where key is k-mer string and value is a list of starting positions
    dna_dict = Dict{String, Vector{Int}}()

    # Populate the dictionary with k-mers from dna1
    # Note: Julia is 1-indexed.
    for i in 1:(length(dna1) - k + 1)
        kmer = dna1[i:i+k-1]
        rev_kmer = reverse_complement(kmer)

        # Add original k-mer
        push!(get!(dna_dict, kmer, Int[]), i)
        # Add reverse complement
        push!(get!(dna_dict, rev_kmer, Int[]), i)
    end

    common_kmers = Set{Tuple{Int, Int}}()

    # Check k-mers in dna2 against the dictionary
    for j in 1:(length(dna2) - k + 1)
        kmer2 = dna2[j:j+k-1]

        # If this k-mer exists in dna1 (or its rev_comp did),
        # add all associated indices from dna1
        if haskey(dna_dict, kmer2)
            for x in dna_dict[kmer2]
                # Storing as (x-1, j-1) if you need 0-indexed output for Stepic
                push!(common_kmers, (x-1, j-1))
            end
        end
    end

    return common_kmers
end

function main()
    # Read the input data
    input_lines = readlines("data/stepic_8d.txt")
    k = parse(Int, strip(input_lines[1]))
    dna1 = strip(input_lines[2])
    dna2 = strip(input_lines[3])

    # Get the shared kmers and sort them for consistency
    common = sort(collect(shared_kmers(k, dna1, dna2)))

    # Format output as "(x, y)" strings
    output_lines = ["($x, $y)" for (x, y) in common]

        # Print result
        println(join(output_lines, "\n"))

        # Save the answer
        mkpath("output")
        write("output/Assignment_08D.txt", join(output_lines, "\n"))
    end

    if abspath(PROGRAM_FILE) == @__FILE__
        main()
    end

