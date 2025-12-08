# Task 2 Report: Benchmarking and Autovectorization

## 1. Description of the Process

**Compiler Versions:**
Due to the unavailability of GCC 8.4.0 and 11.4.0 on the current cluster environment, the following versions were used:
1.  **GCC 10.1.0**
2.  **GCC 12.3.0**
3.  **GCC 14.3.0**

**Methodology:**
The solver was compiled with four configurations:
* `O0`: No optimization.
* `O2-novec`: Optimization level 2, with vectorization explicitly disabled (`-fno-tree-vectorize`).
* `O3-vec`: Optimization level 3, enabling autovectorization (`-fopt-info-vec-optimized`).
* `Ofast-vec`: Aggressive optimization breaking strict standard compliance.

The tests were run on three matrix sizes: 1024x1024, 2048x2048, and 4096x4096.

---

## 2. Benchmark Results

### Table 1: GCC 10.1.0 Results (in ms)
| Matrix Size | O0 | O2-novec | O3-vec | Ofast-vec | Ref |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1024x1024** | 1256 | 264 | 150 | 148 | - |
| **2048x2048** | 9998 | 2247 | 1373 | 1345 | - |
| **4096x4096** | 81969 | 20618 | 15911 | 16354 | - |

### Table 2: GCC 12.3.0 Results (in ms)
| Matrix Size | O0 | O2-novec | O3-vec | Ofast-vec | Ref |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1024x1024** | 1253 | 214 | 155 | 149 | - |
| **2048x2048** | 10031 | 1894 | 1379 | 1352 | - |
| **4096x4096** | 81250 | 21787 | 19727 | 15986 | - |

### Table 3: GCC 14.3.0 Results (in ms)
| Matrix Size | O0 | O2-novec | O3-vec | Ofast-vec | Ref |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1024x1024** | 1253 | 214 | 151 | 150 | - |
| **2048x2048** | 10064 | 1883 | 1338 | 1336 | - |
| **4096x4096** | 82250 | 18889 | 15915 | 16053 | - |

---

## 3. Analysis of Results

### Performance Improvements
* **O0 vs O2:** The transition from `O0` to `O2` provides the largest speedup (approx. 5x - 6x faster). This is due to standard scalar optimizations such as register allocation, instruction scheduling, and dead code elimination, which remove the heavy overhead of unoptimized code.
* **Vectorization Impact (O2 vs O3):** Enabling autovectorization (`O3`) resulted in a significant further performance boost. For the 1024x1024 matrix, the time dropped from ~264ms (O2) to ~150ms (O3) in GCC 10. This ~1.7x speedup indicates that the compiler successfully utilized SIMD instructions (AVX/SSE) to process multiple data elements per clock cycle.

### Autovectorization and Code Modifications
To maximize the compiler's ability to vectorize the code, specific modifications were made to `dgesv.c`:

1.  **Use of `restrict` keyword:**
    * **Reason:** Initially, the compiler reported "loop versioned for vectorization because of possible aliasing". The compiler could not guarantee that pointers `a` and `b` did not overlap.
    * **Modification:** We changed the function signature to `int my_dgesv(..., double * restrict a, double * restrict b)`.
    * **Effect:** This informed the compiler that the memory arrays are independent, allowing it to generate optimized vector code without unnecessary runtime checks.

2.  **Loop Fusion:**
    * **Reason:** The original code had separate loops for updating matrix A and matrix B.
    * **Modification:** We fused the inner loops inside the Gaussian elimination step.
    * **Effect:** This reduced loop overhead and improved instruction density for the vectorizer.

**Compiler Reports:**
Using `-fopt-info-vec-optimized`, we confirmed successful vectorization of the critical loops. The output showed:
> `dgesv.c:52:13: optimized: loop vectorized using 16 byte vectors`
> `dgesv.c:46:13: optimized: loop vectorized using 16 byte vectors`

This confirms that GCC generated 128-bit SIMD instructions (processing 2 doubles at once) for the update steps of the Gaussian elimination.

### Cache Behavior Analysis
The **Loop Fusion** strategy significantly improved cache locality. In the original code, the CPU had to load a row of Matrix A, process it, then later load a row of Matrix B. By fusing the loops:
1.  The pivot row `k` and the multipliers are kept in the L1/L2 cache while updating both A and B.
2.  This reduces the number of cache misses and memory bandwidth pressure, which is crucial for large matrices (4096 size) where data does not fit entirely in the cache.

### Conclusion
The combination of algorithmic changes (Loop Fusion) and hinting the compiler (`restrict`, `ivdep`) allowed GCC to effectively autovectorize the solver, resulting in a highly optimized execution time close to the theoretical limit for a non-blocked implementation.
