#!/usr/bin/env julia

#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
Textbook: Bioinformatics Algorithms: An Active-Learning Approach
Authors: Phillip Compeau & Pavel Pevzner

Problem Title: Protein Translation Problem
=#

# -----------------------------
# RNA Codon Table
# -----------------------------
const RNA_CODON_TABLE = Dict(
    "UUU"=>"F", "UUC"=>"F", "UUA"=>"L", "UUG"=>"L",
    "CUU"=>"L", "CUC"=>"L", "CUA"=>"L", "CUG"=>"L",
    "AUU"=>"I", "AUC"=>"I", "AUA"=>"I", "AUG"=>"M",
    "GUU"=>"V", "GUC"=>"V", "GUA"=>"V", "GUG"=>"V",

    "UCU"=>"S", "UCC"=>"S", "UCA"=>"S", "UCG"=>"S",
    "CCU"=>"P", "CCC"=>"P", "CCA"=>"P", "CCG"=>"P",
    "ACU"=>"T", "ACC"=>"T", "ACA"=>"T", "ACG"=>"T",
    "GCU"=>"A", "GCC"=>"A", "GCA"=>"A", "GCG"=>"A",

    "UAU"=>"Y", "UAC"=>"Y", "UAA"=>"Stop", "UAG"=>"Stop",
    "CAU"=>"H", "CAC"=>"H", "CAA"=>"Q", "CAG"=>"Q",
    "AAU"=>"N", "AAC"=>"N", "AAA"=>"K", "AAG"=>"K",
    "GAU"=>"D", "GAC"=>"D", "GAA"=>"E", "GAG"=>"E",

    "UGU"=>"C", "UGC"=>"C", "UGA"=>"Stop", "UGG"=>"W",
    "CGU"=>"R", "CGC"=>"R", "CGA"=>"R", "CGG"=>"R",
    "AGU"=>"S", "AGC"=>"S", "AGA"=>"R", "AGG"=>"R",
    "GGU"=>"G", "GGC"=>"G", "GGA"=>"G", "GGG"=>"G"
)

# -----------------------------
# Read RNA Sequence
# -----------------------------
input_path = "data/stepic_2a.txt"
rna_seq = open(input_path, "r") do io
    strip(read(io, String))
end

# -----------------------------
# Translate RNA to Protein
# -----------------------------
protein = IOBuffer()

for i in 1:3:length(rna_seq)-2
    codon = rna_seq[i:i+2]
    amino = RNA_CODON_TABLE[codon]
    if amino == "Stop"
        break
    end
    print(protein, amino)
end

s_protein = String(take!(protein))
println(s_protein)
