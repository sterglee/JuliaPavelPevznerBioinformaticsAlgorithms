#!/usr/bin/env julia
"""
Change Problem (Julia version)
"""

# -----------------------------
# Dynamic Programming function
# -----------------------------
function DPChange(amount::Int, coin_list::Vector{Int})
    # Initialize min_coins array: 0 for amount 0, Inf for others
    min_coins = [0; fill(amount ÷ minimum(coin_list) + 1, amount)]

    # Build up the solution using DP
    for m in 1:amount
        for coin in coin_list
            if m >= coin
                min_coins[m+1] = min(min_coins[m+1], min_coins[m-coin+1] + 1)
            end
        end
    end
    return min_coins[amount+1]
end

# -----------------------------
# Main Program
# -----------------------------
filename = "data/stepic_6a.txt"
money, coins = open(filename) do io
    money = parse(Int, readline(io))
    coins = parse.(Int, split(chomp(readline(io)), ","))
    return money, coins
end

# Compute the minimum number of coins
min_number = DPChange(money, coins)

# Print and save the answer
println(min_number)
