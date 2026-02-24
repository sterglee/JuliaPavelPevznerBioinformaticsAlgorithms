"""
Greedy Sorting Algorithm for permutations.
    Returns a list of all intermediate permutations during the sorting process.
    """
    function greedy_sorting(permutation::Vector{Int})
        transformations = Vector{Vector{Int}}()
        perm = copy(permutation)
        n = length(perm)

        # Helper: Reverse and negate a sub-segment from index i to j
        function k_sort!(p, i, j)
            # Reverse the sub-segment and negate each element
            p[i:j] = -reverse(p[i:j])
            return p
        end

        for i in 1:n
            # Goal: Place the value 'i' at index 'i'
            target = i

            # Check if the element at the current position is not what we want
            if perm[i] != target
                # Find where the value 'target' is (ignoring sign)
                # findfirst returns the index where the absolute value matches
                idx = findfirst(x -> abs(x) == target, perm)

                # Step 1: Reverse from current position i to the found index
                k_sort!(perm, i, idx)
                push!(transformations, copy(perm))

                # Step 2: If it's now at the right spot but negative, flip the sign
                # (The Python logic handles this by continuing the loop;
                # we check if perm[i] is now -target)
                if perm[i] == -target
                    k_sort!(perm, i, i)
                    push!(transformations, copy(perm))
                end
            end
        end

        return transformations
    end

    """
    Formats a permutation vector into the string format required by Stepic:
    (+1 -2 -3 +4)
    """
    function format_permutation(perm::Vector{Int})
        formatted = join([x > 0 ? "+$x" : "$x" for x in perm], " ")
            return "($formatted)"
        end

        # --- Main Routine ---

        function main()
            input_path = "data/stepic_8a.txt"
            output_path = "output/Assignment_08A.txt"

            if !isfile(input_path)
                println("Input file not found.")
                return
            end

            # Read input: strip parens and split into integers
            raw_data = read(input_path, String)
            # Regex to extract numbers from something like (+1 -2 -3 +4)
            perm = [parse(Int, m.match) for m in eachmatch(r"[-+]?\d+", raw_data)]

                # Run Greedy Sorting
                reversal_steps = greedy_sorting(perm)

                # Format output
                output_lines = [format_permutation(step) for step in reversal_steps]

                    # Print and Save
                    result_string = join(output_lines, "\n")
                    println(result_string)

                    mkpath("output")
                    write(output_path, result_string)
                end

                main()
