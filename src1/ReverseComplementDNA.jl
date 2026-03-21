#!/usr/bin/env julia
#=
A solution to a programming assignment for the Bioinformatics Algorithms (Part 1) on Coursera.
The associated textbook is Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
=#

include("BioUtils.jl")
using .BioinformaticsUtils: ReverseComplementDNA


dna = strip(read("data/stepic_1b.txt", String))

 println(ReverseComplementDNA(dna))

