#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <chrono>
#include <thread>
#include <atomic>

const int d = 16;
const int NUM_RUNS = 1000;
const int N7 = 16384; 
const int N8 = 65536; 

#define HASH_SIZE (1 << 18)
#define MASK (HASH_SIZE - 1)

struct Bucket {
    int64_t key;
    uint32_t last_run;
    uint32_t c_m2;    
    uint32_t c_clean; 
};

Bucket table[HASH_SIZE] __attribute__((aligned(64)));
int64_t a[d], v[16][4];
int8_t alpha[4] = {-2, -1, 0, 1};
uint32_t current_run = 0;

inline uint32_t hash_func(int64_t val) {
    uint64_t x = (uint64_t)val;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33;
    return (uint32_t)x & MASK;
}

// Global to hold results for the threads
std::atomic<int64_t> global_total_2x(0);

void search_segment(int start_idx, int end_idx, int64_t q) {
    int64_t local_2x = 0;
    int y[8] = {0};
    int64_t s2 = 0;
    int r_m2_count = 0;

    // Fast-forward odometer to start_idx
    for (int i = 0; i < 8; i++) {
        int val_idx = (start_idx >> (2 * i)) & 3;
        y[i] = val_idx;
        s2 += v[i + 8][y[i]];
        if (alpha[y[i]] == -2) r_m2_count++;
    }

    for (int i = start_idx; i < end_idx; i++) {
        int64_t m2 = s2 % q; if (m2 < 0) m2 += q;
        bool r_has_m2 = (r_m2_count > 0);

        for (int k = 0; k < 4; k++) {
            int8_t x0 = alpha[k];
            int64_t target = x0 - m2;
            if (target < 0) target += q; else if (target >= q) target -= q;

            uint32_t h = hash_func(target);
            while (table[h].last_run == current_run) {
                if (table[h].key == target) {
                    if (x0 == -2 || r_has_m2) local_2x += 2 * (table[h].c_m2 + table[h].c_clean);
                    else local_2x += 2 * table[h].c_m2 + table[h].c_clean;
                }
                h = (h + 1) & MASK;
            }
        }

        // Odometer update
        for (int j = 0; j < 8; j++) {
            int idx = j + 8;
            if (++y[j] < 4) {
                if (alpha[y[j] - 1] == -2) r_m2_count--;
                if (alpha[y[j]] == -2) r_m2_count++;
                s2 += (v[idx][y[j]] - v[idx][y[j] - 1]);
                break;
            } else {
                if (alpha[3] == -2) r_m2_count--;
                if (alpha[0] == -2) r_m2_count++;
                s2 += (v[idx][0] - v[idx][3]);
                y[j] = 0;
            }
        }
    }
    global_total_2x += local_2x;
}

int solve() {
    current_run++;
    int64_t q = a[0];
    global_total_2x = 0;

    // 1. Left Half (Single Threaded - Write heavy)
    int x[7] = {0};
    int64_t s1 = 0;
    int l_m2_count = 0;
    for(int i=1; i<=7; i++) {
        s1 += v[i][0];
        if (alpha[0] == -2) l_m2_count++;
    }

    for(int i=0; i<N7; i++) {
        int64_t m1 = s1 % q; if (m1 < 0) m1 += q;
        bool has_m2 = (l_m2_count > 0);
        uint32_t h = hash_func(m1);
        
        static std::atomic<int64_t> collisions{0};
        // ... inside Phase 1 loop, after hash_func(m1) ...
        while (table[h].last_run == current_run) {
            if (table[h].key == m1) { 
                if (++collisions < 10) printf("Duplicate sum %ld at left_idx %d\n", m1, i);
                goto next_i; 
            }
            h = (h + 1) & MASK;
        }

        while (table[h].last_run == current_run) {
            if (table[h].key == m1) {
                if (has_m2) table[h].c_m2++; else table[h].c_clean++;
                goto next_i;
            }
            h = (h + 1) & MASK;
        }
        table[h].key = m1; table[h].last_run = current_run;
        table[h].c_m2 = has_m2 ? 1 : 0;
        table[h].c_clean = has_m2 ? 0 : 1;
    next_i:
        for(int j=0; j<7; j++) {
            int idx = j + 1;
            if(++x[j] < 4) {
                if (alpha[x[j]-1] == -2) l_m2_count--;
                if (alpha[x[j]] == -2) l_m2_count++;
                s1 += (v[idx][x[j]] - v[idx][x[j]-1]);
                break;
            } else {
                if (alpha[3] == -2) l_m2_count--;
                if (alpha[0] == -2) l_m2_count++;
                s1 += (v[idx][0] - v[idx][3]);
                x[j] = 0;
            }
        }
    }

    // 2. Right Half (Split into 2 threads)
    std::thread t1(search_segment, 0, N8/2, q);
    search_segment(N8/2, N8, q); // Main thread does the other half
    t1.join();

    return (int)(global_total_2x / 2);
}

int main(int argc, char** argv) {
    if (argc < 18) return 1;
    for (int i = 0; i < d; i++) a[i] = atoll(argv[i+2]);
    for (int i = 0; i < 16; i++) {
        for(int k=0; k<4; k++) v[i][k] = (int64_t)alpha[k] * a[i];
    }
    auto start = std::chrono::high_resolution_clock::now();
    int res = 0;
    for (int i = 0; i < NUM_RUNS; i++) res = solve();
    auto end = std::chrono::high_resolution_clock::now();
    printf("Benchmark: %.2f us | Unique: %d\n", (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / NUM_RUNS, res);
    return 0;
}

