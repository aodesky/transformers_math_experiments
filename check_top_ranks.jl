#!/usr/bin/env julia

"""
Check ranks of top N rational values from search results
Uses GP/PARI (fast) - NO SageMath
"""

using Printf

include("problem_mestre_rank12.jl")

function check_top_results(search_output_file::String, top_n::Int=50)
    """
    Read top N results from search_output file and check their ranks
    """
    println("="^80)
    println("Checking ranks of top $top_n results from $search_output_file")
    println("Using GP/PARI (fast, not Sage)")
    println("="^80)
    println()

    # Read top N lines
    rationals = String[]
    open(search_output_file, "r") do f
        for i in 1:top_n
            line = readline(f)
            if !isempty(line)
                push!(rationals, strip(line))
            else
                break
            end
        end
    end

    println("Found $(length(rationals)) rationals to check\n")

    # Check each one
    results = []
    for (i, obj) in enumerate(rationals)
        print("[$i/$(length(rationals))] Checking t=$obj ... ")
        flush(stdout)

        try
            num, den = parse_rational(obj)

            # Check discriminant
            disc = eval_discriminant(num, den)
            if abs(disc) < 1e-10
                println("❌ Zero discriminant")
                continue
            end

            # Compute rank (with 30s timeout)
            rank_lower = -1
            a4_val = eval_a4(num, den)
            a6_val = eval_a6(num, den)

            a4_num = numerator(a4_val)
            a4_den = denominator(a4_val)
            a6_num = numerator(a6_val)
            a6_den = denominator(a6_val)

            gp_script = """
            a4 = $a4_num/$a4_den;
            a6 = $a6_num/$a6_den;
            E = ellinit([0, 0, 0, a4, a6]);
            if (E == 0, print("INVALID"), r = ellrank(E); print(r[1]));
            quit();
            """

            temp_file = tempname() * ".gp"
            open(temp_file, "w") do f
                write(f, gp_script)
            end

            result_str = read(pipeline(`timeout 30 gp -q -s 200000000 $temp_file`, stderr=devnull), String)
            rm(temp_file)
            result_str = strip(result_str)

            if result_str != "INVALID" && !isempty(result_str)
                rank_lower = parse(Int, result_str)
            end

            # Compute conductor
            conductor = compute_conductor_mestre(num, den, check_rank=false)

            if rank_lower >= 12
                println("✓ Rank ≥ $rank_lower, Conductor = $conductor")
                push!(results, (obj, rank_lower, conductor, true))
            elseif rank_lower >= 0
                println("⚠️  Rank = $rank_lower (< 12), Conductor = $conductor")
                push!(results, (obj, rank_lower, conductor, false))
            else
                println("❌ Rank computation failed")
            end

        catch e
            println("❌ Error: $e")
        end
    end

    # Summary
    println()
    println("="^80)
    println("SUMMARY")
    println("="^80)

    rank_12_plus = filter(r -> r[4], results)
    println("Found $(length(rank_12_plus)) curves with rank ≥ 12:")
    println()

    if !isempty(rank_12_plus)
        println(@sprintf("%-15s | %-10s | %s", "t", "Rank", "Conductor"))
        println("-"^80)
        for (t, r, c, _) in rank_12_plus
            println(@sprintf("%-15s | %-10d | %s", t, r, c))
        end

        println()
        println("Best conductor among rank ≥ 12 curves:")
        best = minimum(r[3] for r in rank_12_plus)
        best_t = filter(r -> r[3] == best, rank_12_plus)[1][1]
        println("  t = $best_t")
        println("  Conductor = $best")
    else
        println("⚠️  No curves with rank ≥ 12 found in top $top_n results")
    end
    println("="^80)
end

# Main
if length(ARGS) < 1
    println("Usage: julia check_top_ranks.jl <search_output_file> [top_n]")
    println("Example: julia check_top_ranks.jl checkpoint/debug/xyz/search_output_1.txt 100")
    exit(1)
end

search_file = ARGS[1]
top_n = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50

check_top_results(search_file, top_n)
