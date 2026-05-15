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

## Tools

### 1. makesievebase
Before sieving, this utility generates the required factor base and roots 
for the polynomial. It is parallelized with OpenMP for rapid generation.

**Usage:**
./makesievebase [polyfile] [pmax] [output_sb] [threads]
```bash
./makesievebase rsa1024b.poly 2147483647 rsa1024b.M31.sb 128
```
## 2. hdlcs
The primary sieving binary. It performs lattice enumeration, linear sweeping
of the 4GB arrays, and cofactorization of potential relations.

### Parameters:
* inputpoly: Polynomial file (N/skew/C/Y format).
* sievebase: Pre-computed sieve base from makesievebase.
* d: Sieving dimension (fixed at 16).
* Amax/Bmax: Bounds for the ideal generator coefficients.
* N: Number of workunits to process.
* pmin/pmax: Range of sieving primes.
* th0/th1: Log-sum thresholds for Algebraic and Rational sides.
* lpb: Large Prime Bound.ecmpbits: Bit-length for ECM (typically 11).
* bb: Bits in lattice coefficient range (typically 2).
* seed: Random seed for reproducibility.

### Example Execution:
```bash
$ time taskset -c 0,1,2,3 ./hdlcs rsa1024b.poly rsa1024b.M31.sb 16 5000 100000000 1 2000 100000000 90 70 549755813887 11 2 12345 | tee test001.rels
```

## Implementation Details
hdlcs.cc The code is organized to minimize synchronization overhead and maximize
memory throughput.
### MITM_Workspace:
Encapsulates the 4GB sieve array, bucket-sieving buffers, and the
shared hash table used for enumeration.
### 4-Thread Architecture:
* Threads 0 & 1: Algebraic side MITM.
* Threads 2 & 3: Rational side MITM.
### Enumeration Phases:
* Phase 1: One thread populates a hash table with the "Left Half"of the 16D orthotope.
* Phase 2: Both threads search the "Right Half" in parallel,emitting hits to
region-specific buckets.
### Linear Sweep:
A multi-threaded scan of the 4GB arrays to identify indiceswhere both
sides exceed the specified thresholds.
### Cofactorization:
Uses Pollard's P-1 and EECM (Edwards Elliptic Curve Method) to find
actual relations from potential hits.
### Performance & Projections
The implementation achieves stable performance across both low and high prime ranges.
Full Range Projection: For a search up to $2^{31}-1$ (approx. 105 million primes),
the projected runtime is 12-13 hours on a 4-core allocation.
### Efficiency:
The code utilizes a memset clear of the 4GB arrays onlyonce per workunit,
ensuring that the bulk of the time is spent on high-value MITM calculations rather
than memory management.

## Producing relations
The main binary hdlcs is now capable of finding relations for RSA-1024.
Using the polynomial rsa1024b.poly in this repository, some test sieving was carried out.
The run time is not good but there is a lot of optimization left to do.

Here is the output of one instance which found a relation:

```bash
# ./hdlcs rsa1024b.poly rsa1024b.M31.sb 16 5000 100000000 1 1000 200000000 100 65 4398046511103 11 2 1812476
# Loading sieve base data...done.
# Computing small prime array...done.
# ========== Unit 1 / 1 ==========
# 3210*x + 20658519
# 3208*x + 5297027
# 2816*x + 78117791
# 1800*x + 89697394
# 1067*x + 10483530
# 3504*x + 77052357
# 2821*x + 98041160
# 2706*x + 20292919
# 1306*x + 15952415
# 1683*x + 46678970
# 2141*x + 65469526
# 4083*x + 30644137
# 800*x + 96797886
# 901*x + 18797244
# Clearing 2x sieve arrays (4GB each)...done.
# Launching 4-thread execution for HDLCS parallel MITM...
# HDLCS Enumeration Finished! Time: 13258.8s
# Executing concurrent multi-threaded linear sweep...
# Linear Sweep Finished! Time: 0.845064s
# 31 potential relations found.
# Starting cofactorization...
1172,-25655295:2,2,2,2,2,3,5,7,7,d,2b,ef,2190b,244469,6b69d7,8895ad,21b5345,3292205,6929beb,3713fe45,604adb51db:b,1bfb5,2225f,3a5577,422b37,7588ac9,a713f8f,e68363e23
# Cofactorization took 48.2435s
# 1 actual relations found.
```

A typical slurm script running on a 128 core node might look like this:

```bash
#!/bin/bash -l
#SBATCH -N 1
#SBATCH --time=05:00:00
#SBATCH --account=NNNNNNN
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH --cpus-per-task=128
#SBATCH --mail-type=END
#SBATCH --mail-user=name@domain.com

cd $SLURM_SUBMIT_DIR

module load GCC
module load GMP
module load GDB

# Array of 32 seeds
seeds=(5761819 2247845 6872380 9702806 6256416 6243731 8783914 6506512 \
       8051602 6362647 6898766 4927732 9968914 3487059 8864868 5026344 \
       4648617 2414845 1445631 7979571 4858867 138917 2797155 5541244 \
       3709493 7075950 107202 9433127 9009015 8816339 1812476 5157733)

# Parameters
POLY="rsa1024b.poly"
SB="rsa1024b.M31.sb"

# Loop through 0-31 to launch 32 instances
for i in {0..31}; do
    SEED=${seeds[$i]}

    # Calculate CPU affinity range (e.g., 0-3, 4-7, 8-11...)
    CPU_START=$((i * 4))
    CPU_END=$((CPU_START + 3))
    CPUS="$CPU_START-$CPU_END"

    echo "Launching instance $i with seed $SEED on CPUs $CPUS"

    # Run in background (&). Each output is logged to a unique file.
    taskset -c "$CPUS" ./hdlcs "$POLY" "$SB" 16 5000 100000000 1 1000 \
    200000000 100 65 4398046511103 11 2 "$SEED" > "run_${SEED}.rels" 2>&1 &
done

# Wait for all background processes to finish
wait

echo "All 32 instances completed."
```

### fast39.cc
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

## Example
```bash
$ time taskset -c 0,1 ./fast39 2 97668217 27746440 9402514 8743350 55200183 76204135 25242614 58741714 77664091 71711529 55743037 60840359 117016 16068556 93981235 18969686
Benchmark: 1003.15 us | Unique: 53

real    0m1.009s
user    0m1.735s
sys     0m0.043s

$ time taskset -c 0,1 ./fast39 2 9913 8092 5930 8439 3512 6846 2196 5202 3341 3258 9673 3544 5953 7008 8399 7757
Benchmark: 1602.79 us | Unique: 430590

real    0m1.607s
user    0m2.921s
sys     0m0.057s
```

## Disclaimer

This is experimental research code. It is not production-hardened and may
require adaptation for large-scale NFS jobs.

---

## License

(Todo)

---

## Author

Oisin Robinson

