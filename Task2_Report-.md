

Environment:

The benchmarks were executed on a local Windows machine using GCC 15.2.0 (MinGW-w64) because the university cluster and BLAS libraries were not available.
The code and results are stored in the GitHub repository.


Benchmark Results:

|        Size   | O0          | O2-novec    | O3-vec      | Ofast-vec   |
|               |             |             |             |             |
| **1024x1024** | 969 ms      | 287 ms      | **226 ms**  | 253 ms      |
| **2048x2048** | 7873 ms     | 3166 ms     | **2511 ms** | 2504 ms     |
| **4096x4096** | 61495 ms    | 24957 ms    | **21485 ms**| 23085 ms    |

Performance Analysis

O0 → O2: This gives the biggest speedup.
The compiler removes unnecessary operations and uses registers better.

O2 → O3: Enabling vectorization gives around 10–15% improvement, especially for the largest matrix.

O3 → Ofast: -Ofast did not improve performance on this machine.
-O3 -march=native was the best option.

Cache Behaviour

The execution time grows almost exactly like O(N³).
When the matrix size doubles, the time becomes 8× larger.

This means the program has good cache locality, because the code accesses rows in order (row-major), which fits well in cache.


Vectorization Efficiency

The comparison between `O2-novec` (No Vectorization) and `O3-vec` (Vectorization Enabled) highlights the success of the code modifications.

- For the largest matrix (4096), the time improved from **24.9s** (O2) to **21.5s** (O3).
- This **14% speedup** indicates that GCC successfully utilized SIMD instructions (AVX) to process multiple double-precision elements per clock cycle.

Autovectorization

The main inner loop (updating row elements) was vectorized 

Some loops with dependencies (pivot search) were not vectorized 

To help vectorization, I made the loop simpler and added restrict to pointers.
This reduced execution time (e.g., 24.9s → 21.5s on the 4096 matrix)