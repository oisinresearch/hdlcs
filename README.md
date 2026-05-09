# hdlcs: High Dimensional Linear Combination Sieve

A high-performance implementation of a High Dimensional Linear Combination 
Sieve (HDLCS) designed for the Number Field Sieve (NFS).

## Overview

The `hdlcs` project introduces a novel sieving strategy for the Number 
Field Sieve (NFS) that utilizes linear combinations of a set of linear 
polynomials rather than traditional lattice sieving directly over (a,b)
for polynomials a+bx. By operating in a fixed high dimension (d=16), the sieve 
targets a compact search space that fits within 4GB of memory.

This approach offers a theoretical 10x to 20x throughput improvement over 
traditional GNFS sieving methods (e.g., I=17 or I=18), which typically 
require 40GB to 80GB per siever. While hardware constraints (RAM per core) 
may necessitate a reduction in core count, the performance is reclaimed 
using a dual-threaded Meet-in-the-Middle (MITM) enumeration strategy.

## Technical Concept

The project produces a "sieving ideal" by finding linear combinations of a 
vector of linear polynomials with skewed sizes. We ensure divisibility by 
all factor base primes p by sieving for a lattice modulo each p.

### Lattice Configuration
Each sieving prime p defines a knapsack-style lattice. Given a root r of 
the NFS polynomial f mod p, and a vector of linear polynomials 
{Li(x) = ai*x + bi}, the lattice is constructed as follows:

* The prime p is placed in the top-left entry.
* The value r is used as a root.
* The remaining first-row vectors are defined by ai*r + bi mod p.

While traditional approaches might use LLL (Lenstra-Lenstra-Lovász) to 
reduce this lattice, LLL is computationally prohibitive in 16 dimensions. 
Instead, `hdlcs` employs a highly optimized MITM enumeration to discover
(all) vectors in this lattice that are contained in the [-2, 1]^16 orthotope.
Note that this is not a simple translation to 16 dimensional sieving ideals
for which the norms would be disastrous.  The high dimension space is only
used to discover a set of weights for the linear polynomials - the weight sum
then gives the ultimate sieving ideal and we have arranged the norm to be
divisible by many primes.

### The Orthotope and Memory Efficiency
We fix the dimension at d=16 with coefficients in the orthotope [-2, 1]^16.

* **Space Complexity:** This search space occupies exactly 4GB, allowing 
  an entire siever to reside in the RAM typically allocated to a single 
  compute core on modern HPC clusters like Meluxina.
* **Hardware Alignment:** To balance the memory-to-core ratio, we utilize 
  2 threads per siever. This configuration speeds up the MITM enumeration 
  and maximizes throughput on nodes where memory is a bottleneck.

## Performance & Feasibility

The `fast39.cc` utility serves as a feasibility test to determine if it is 
possible to enumerate approximately Pi(2^31) sieving prime lattices into a 
4GB sieve array. This involves processing roughly 100 million lattices.

* **Target Speed:** Ideally < 1ms per lattice enumeration.
* **Results:** The current implementation achieves near 1ms performance 
  per enumeration. This speed remains consistent for both low vector 
  counts (< 100, 31 bit prime) and high vector counts (> 400k, 12 bit prime).

## Implementation Details

`fast39.cc` is the result of intensive optimization:

* **Dual-Threaded MITM:** Parallelizes the heaviest part of the 
  enumeration (the Right-Half search).
* **Incremental Odometer:** Updates running sums incrementally to 
  minimize expensive modular arithmetic.
* **Symmetry Exploitation:** Categorizes vectors (clean vs. containing -2) 
  to exploit mathematical symmetry, effectively halving the search space 
  without losing correctness.

## Usage

### Compilation
We use pthreads:
```bash
g++ -O3 -march=native -flto=auto -pthread -falign-functions=32 \
    -fno-plt fast39.cc -o fast39
```

## Execution
The invocation format is:

```bash
taskset -c 0,1 ./fast39 [B] [q] [a2] [a3] ... [a16]
```

B: The orthotope bound. For the [-2, 1]^16 space, use 2.

q: The modulus (sieving prime).

a2...a16: The lattice basis elements.

## Disclaimer

This is experimental research code. It is not production-hardened and may
require adaptation for large-scale NFS jobs.

---

## License

(Todo)

---

## Author

Oisin Robinson

