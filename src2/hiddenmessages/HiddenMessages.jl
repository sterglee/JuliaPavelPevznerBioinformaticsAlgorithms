# Bioinformatics Algorithms in Julia
# Copyright (C) 2016  Jason Mar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

module HiddenMessages

export patternCount, frequentWords, reverseComplement

"""
patternCount(text::String, pattern::String) -> Int

Count the number of times `pattern` appears in `text`.
"""
function patternCount(text::String, pattern::String)::Int
    count = 0
    plen = length(pattern)
    for i in 1:(length(text) - plen + 1)
        if text[i:i+plen-1] == pattern
            count += 1
        end
    end
    return count
end

"""
frequentWords(text::String, k::Int) -> Vector{String}

Return all most frequent k-mers in `text`.
"""
function frequentWords(text::String, k::Int)::Vector{String}
    kmerCounts = kmerFrequency(text, k)
    maxCount = maximum(values(kmerCounts))
    most_frequent = [kmer for (kmer, count) in kmerCounts if count == maxCount]
        sort(most_frequent)
    end

    # Helper function: returns a dictionary mapping k-mers to their counts
    function kmerFrequency(text::String, k::Int)::Dict{String, Int}
        kmers = getKmers(text, k)
        counts = Dict{String, Int}()
        for kmer in kmers
            counts[kmer] = get(counts, kmer, 0) + 1
        end
        return counts
    end

    # Helper function: returns all k-mers in text (including duplicates)
    function getKmers(text::String, k::Int)::Vector{String}
        kmers = Vector{String}(undef, length(text) - k + 1)
        for i in 1:(length(text) - k + 1)
            kmers[i] = text[i:i+k-1]
        end
        return kmers
    end

    # Complement a nucleotide character
    function complement(p::Char)::Char
        c = uppercase(p)
        if c == 'A'
            return 'T'
            elseif c == 'C'
            return 'G'
            elseif c == 'G'
            return 'C'
            elseif c == 'T'
            return 'A'
        else
            return 'N'
        end
    end

    """
    reverseComplement(pattern::String) -> String

    Return the reverse complement of a DNA string `pattern`.
    """
    function reverseComplement(pattern::String)::String
        chars = [complement(c) for c in reverse(pattern)]
            return join(chars)
        end

    end # module
