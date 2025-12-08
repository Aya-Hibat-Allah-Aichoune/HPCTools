#!/bin/bash

# Sizes required
SIZES="1024 2048 4096"

# Check version
echo "----------------------------------------------------"
echo "Running Benchmark with GCC version: $(gcc --version | head -n 1)"
echo "----------------------------------------------------"
echo "Size      | Flag           | Time (ms)"
echo "----------|----------------|----------"

# Function to compile and run
run_test() {
    FLAG_NAME=$1
    FLAGS=$2
    
    # Compile (Removed /dev/null to see errors if any)
    gcc $FLAGS main.c dgesv.c timer.c -o bench_exec
    
    # Check if compilation succeeded
    if [ ! -f ./bench_exec ]; then
        echo "Error: Compilation failed for $FLAG_NAME"
        return
    fi

    # Run for each size
    for N in $SIZES; do
        # Grep output: "Time taken by my_dgesv: 151 ms"
        # We need field $5 (the number)
        OUTPUT=$(./bench_exec $N | grep "Time" | awk '{print $5}')
        echo "$N      | $FLAG_NAME             | $OUTPUT"
    done
}

# 1. O0 (No optimization)
run_test "O0" "-O0"

# 2. O2-novec (Optimization Level 2, No Vectorization)
run_test "O2-novec" "-O2 -fno-tree-vectorize"

# 3. O3-vec (Optimization Level 3, Auto-Vectorization enabled)
run_test "O3-vec" "-O3 -fopt-info-vec-optimized"

# 4. Ofast-vec (Aggressive Optimization)
run_test "Ofast-vec" "-Ofast -fopt-info-vec-optimized"

# Clean up
rm -f bench_exec
