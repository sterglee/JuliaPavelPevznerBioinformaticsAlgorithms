# -----------------------------
# Amino Acid Integer Mass Table
# -----------------------------
const AMINO_ACID_MASS = Dict(
    'G'=>57,  'A'=>71,  'S'=>87,  'P'=>97,
    'V'=>99,  'T'=>101, 'C'=>103,'I'=>113,
    'L'=>113, 'N'=>114,'D'=>115,'K'=>128,
    'Q'=>128, 'E'=>129,'M'=>131,'H'=>137,
    'F'=>147, 'R'=>156,'Y'=>163,'W'=>186
)

# -----------------------------
# Generate Cyclic Spectrum

function cyclospectrum(peptide::AbstractString, AMINO_ACID_MASS::Dict{Char, Int})
    n = length(peptide)
    peptide_mass = sum(AMINO_ACID_MASS[c] for c in peptide)
    spectrum = [0, peptide_mass]

    # Repeat the peptide to handle cyclic subpeptides
    double_peptide = repeat(peptide, 2)

    for sub_len in 1:(n-1)
        for start in 1:n
            subpep = SubString(double_peptide, start, start + sub_len - 1)
            subpep_mass = sum(AMINO_ACID_MASS[c] for c in subpep)
            push!(spectrum, subpep_mass)
        end
    end

    return sort(spectrum)
end


# -----------------------------
# Main
# -----------------------------
peptide = "RRNQKRGCLSQQCFL"
cyclospec = cyclospectrum(peptide, AMINO_ACID_MASS)

# Convert to strings and join for output
result = join(cyclospec, " ")
println(result)
