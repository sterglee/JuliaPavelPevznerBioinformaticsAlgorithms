using Printf, LinearAlgebra

function ProteinWeightDict()
    return Dict(
        'G' => 57, 'A' => 71, 'S' => 87, 'P' => 97,
        'V' => 99, 'T' => 101, 'C' => 103, 'I' => 113,
        'L' => 113, 'N' => 114, 'D' => 115, 'K' => 128,
        'Q' => 128, 'E' => 129, 'M' => 131, 'H' => 137,
        'F' => 147, 'R' => 156, 'Y' => 163, 'W' => 186
    )
end

function append_char(add_list::Vector{String}, add_chars::Vector{String})
    newlist = String[]
    for item in add_list
        for ch in unique(add_chars)
            push!(newlist, item * ch)
        end
    end
    return newlist
end

function spectrum(peptide::String, weight::Dict{Char, Int})
    n = length(peptide)
    spec = [0, sum(weight[c] for c in peptide)]
    for i in 1:(n-1)
        for j in 1:(n-i+1)
            push!(spec, sum(weight[c] for c in peptide[j:j+i-1]))
        end
    end
    sort!(spec)
    return spec
end

# -----------------------------
# Main
# -----------------------------
cyclospec = parse.(Int, split("0 97 113 113 131 131 137 147 156 163 226 234 244 253 268 269 278 294 310 357 365 366 382 390 407 415 441 441 479 503 512 513 520 521 554 572 578 610 616 634 667 668 675 676 685 709 747 747 773 781 798 806 822 823 831 878 894 910 919 920 935 944 954 962 1025 1032 1041 1051 1057 1057 1075 1075 1091 1188"))

weight = ProteinWeightDict()
L = length(cyclospec)
n = round(Int, (sqrt(4*L - 7) + 1) / 2)

# Extract first n proteins
protein = Int[]
i = 2
weight_values = collect(values(weight))
while length(protein) < n && i <= L
    if cyclospec[i] in weight_values
        push!(protein, cyclospec[i])
    end
    i += 1
end

# Map weights back to protein letters
names = String[]
for w in protein
    key = first(k for (k,v) in weight if v == w)
    push!(names, string(key))
end

# Generate sequences
seq = append_char(names, names)
for repeat in 1:(n-1)
    seq = [s for s in unique(seq) if Set(spectrum(s, weight)) ⊆ Set(cyclospec)]
    if repeat != n-1
        seq = append_char(seq, names)
    end
end

# Convert to weight format
cyclopeptide_sequence = [join([string(weight[Char(c)]) for c in s], "-") for s in seq]

# Output
println(join(cyclopeptide_sequence, " "))

