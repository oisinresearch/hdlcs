/*
 * HDLCS: High Dimensional Linear Combination Sieve
 * * An optimized MITM lattice enumeration approach utilizing fixed 16D orthotopes.
 * Designed for strict 4-thread execution (2 per side) and cache-friendly bucket sieving.
 */

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <gmpxx.h>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stack>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <chrono>
#include <memory>

#include "intpoly.h"
#include "mpz_poly.h"

using namespace std;

#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 uint128_t;
typedef __int128 int128_t;
#else
#error "128-bit integers are required for this NFS build."
#endif

int128_t MASK64;

// ==================== CONSTANTS ====================
const int SIEVE_D = 16;
const uint64_t ARRAY_SIZE = (1ULL << 32);     // 4GB array space
const int NUM_REGIONS = 256;                  // For bucket sieving
const int REGION_SHIFT = 24;                  // 32-bit idx >> 24 yields 0-255 region
const int BUCKET_CAPACITY = 1024;
const uint32_t HASH_SIZE = (1 << 19);
const uint32_t HASH_MASK = (HASH_SIZE - 1);

// ==================== DATA STRUCTURES ====================
struct SieveSide {
    int k;
    vector<int> p;
    vector<int> r;
    vector<int> n;
    vector<int> r_offset;
    uint8_t threshold;
};

// Simple Thread Barrier
struct SpinBarrier {
    std::atomic<int> count{0};
    std::atomic<int> generation{0};
    int num_threads;
    
    SpinBarrier(int n) : num_threads(n) {}
    
    void wait() {
        if (num_threads <= 1) return; // Immediate return for debug/single-thread
        int gen = generation.load(std::memory_order_acquire);
        if (count.fetch_add(1, std::memory_order_acq_rel) == num_threads - 1) {
            count.store(0, std::memory_order_release);
            generation.fetch_add(1, std::memory_order_release);
        } else {
            while (generation.load(std::memory_order_acquire) == gen) {
                std::this_thread::yield();
            }
        }
    }
};

struct LeftNode {
    uint32_t coords;
    uint32_t next;
};

// OPTIMIZATION: logp is removed; entries strictly hold memory locations.
struct Bucket {
    uint32_t idx[BUCKET_CAPACITY];
    int count = 0;
};

struct MITM_Workspace {
    uint8_t* array;
    std::mutex region_locks[NUM_REGIONS];
    
    // Hash Table for Left Half (Shared between thread 0 and 1)
    uint32_t head[HASH_SIZE];
    uint32_t last_run[HASH_SIZE];
    int64_t keys[HASH_SIZE];
    
    LeftNode nodes[16384 + 1];
    uint32_t node_count = 1;
    uint32_t current_run = 0;
    
    SpinBarrier barrier;
    
    // Thread-local buckets [thread_id][region]
    Bucket buckets[2][NUM_REGIONS];
    
#ifdef DEBUG
    MITM_Workspace() : barrier(1) { // Single thread barrier for debug
#else
    MITM_Workspace() : barrier(2) {
#endif
        array = (uint8_t*)calloc(ARRAY_SIZE, 1);
        if (!array) {
            cerr << "CRITICAL ERROR: Failed to allocate 4GB sieve array." << endl;
            exit(1);
        }
        memset(last_run, 0, sizeof(last_run));
    }
    
    ~MITM_Workspace() {
        free(array);
    }
    
    // OPTIMIZATION: Flush a single bucket locally when full, saturating values.
    void emit_to_bucket(int thread_id, uint32_t full_idx, uint8_t logp_val) {
        int region = full_idx >> REGION_SHIFT;
        Bucket& b = buckets[thread_id][region];
        b.idx[b.count] = full_idx;
        b.count++;
        
        if (b.count == BUCKET_CAPACITY) {
            std::lock_guard<std::mutex> lock(region_locks[region]);
            for (int i = 0; i < BUCKET_CAPACITY; i++) {
                uint32_t target_idx = b.idx[i];
                if (array[target_idx] < 255) {
                    uint16_t new_val = (uint16_t)array[target_idx] + logp_val;
                    array[target_idx] = (new_val > 255) ? 255 : (uint8_t)new_val;
                }
            }
            b.count = 0;
        }
    }
    
    // Used when logp changes or when shutting down the side
    void flush_all_buckets(int thread_id, uint8_t logp_val) {
        for (int region = 0; region < NUM_REGIONS; region++) {
            Bucket& b = buckets[thread_id][region];
            if (b.count > 0) {
                std::lock_guard<std::mutex> lock(region_locks[region]);
                for (int i = 0; i < b.count; i++) {
                    uint32_t target_idx = b.idx[i];
                    if (array[target_idx] < 255) {
                        uint16_t new_val = (uint16_t)array[target_idx] + logp_val;
                        array[target_idx] = (new_val > 255) ? 255 : (uint8_t)new_val;
                    }
                }
                b.count = 0;
            }
        }
    }
};

inline uint32_t hash_func(int64_t val) {
    uint64_t x = (uint64_t)val;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33;
    return (uint32_t)x & HASH_MASK;
}

// ==================== FORWARD DECLARATIONS ====================
int64_t rel2A(int d, mpz_ptr* Ai, int64_t reli, int bb);
int64_t rel2B(int d, mpz_ptr* Bi, int64_t reli, int bb);
inline int64_t gcd(int64_t a, int64_t b);
void GetlcmScalar(int B, mpz_t S, int* primes, int nump);
bool PollardPm1(mpz_ptr N, mpz_t S, mpz_t factor);
bool PollardPm1_mpz(mpz_t N, mpz_t S, mpz_t factor);
bool PollardPm1_int128(int128_t N, mpz_t S, int128_t &factor, int64_t, int64_t);
bool EECM(mpz_ptr N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0);
bool EECM_mpz(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0);
bool EECM_int128(int128_t N, mpz_t S, int128_t &factor, int d, int a, int X0, int Y0, int Z0, int64_t, int64_t);

// ==================== I/O HELPERS ====================

void parse_polynomial(const char* filename, mpz_poly f0, mpz_poly f1, double &skew, int &degf, int &degg) {
    ifstream file(filename);
    if (!file.is_open()) { cerr << "Error: Could not open poly file." << endl; exit(1); }
    string line;
    degf = 0; degg = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        size_t colon = line.find(':');
        if (colon == string::npos) continue;
        string key = line.substr(0, colon);
        string val_str = line.substr(colon + 1);
        val_str.erase(0, val_str.find_first_not_of(" \t"));
        if (key == "skew") {
            skew = stod(val_str);
        } else if (key[0] == 'c') {
            int idx = stoi(key.substr(1));
            degf = max(degf, idx);
            mpz_t v; mpz_init_set_str(v, val_str.c_str(), 10);
            mpz_poly_setcoeff(f0, idx, v);
            mpz_clear(v);
        } else if (key[0] == 'Y') {
            int idx = stoi(key.substr(1));
            degg = max(degg, idx);
            mpz_t v; mpz_init_set_str(v, val_str.c_str(), 10);
            mpz_poly_setcoeff(f1, idx, v);
            mpz_clear(v);
        }
    }
    f0->deg = degf; f1->deg = degg;
}

void load_factor_base(const char* filename, SieveSide* sides, uint8_t* th) {
    ifstream fbfile(filename);
    if (!fbfile.is_open()) { cerr << "Error: Could not open FB file." << endl; exit(1); }
    string line;
    getline(fbfile, line); 
    for (int s = 0; s < 2; s++) {
        getline(fbfile, line);
        sides[s].k = stoi(line);
        sides[s].threshold = th[s];
        int offset = 0;
        for (int i = 0; i < sides[s].k; i++) {
            getline(fbfile, line);
            stringstream ss(line);
            string val;
            getline(ss, val, ',');
            int p = stoi(val);
            sides[s].p.push_back(p);
            int count = 0;
            while (getline(ss, val, ',')) {
                sides[s].r.push_back(stoi(val));
                count++;
            }
            sides[s].n.push_back(count);
            sides[s].r_offset.push_back(offset);
            offset += count;
        }
    }
}

// ==================== RELATIONS & COFACTORING ====================

int64_t rel2A(int d, mpz_ptr* Ai, int64_t reli, int bb) {
    int64_t BB = 1LL << bb, hB = 1LL << (bb - 1), A = 0;
    for (int j = 0; j < d - 2; j++) A += mpz_get_si(Ai[j]) * (((reli >> (bb * j)) & (BB - 1)) - hB);
    return A + (((reli >> (bb * (d - 2))) & (BB - 1)) - hB);
}

int64_t rel2B(int d, mpz_ptr* Bi, int64_t reli, int bb) {
    int64_t BB = 1LL << bb, hB = 1LL << (bb - 1), B = 0;
    for (int j = 0; j < d - 2; j++) B += mpz_get_si(Bi[j]) * (((reli >> (bb * j)) & (BB - 1)) - hB);
    return B - (((reli >> (bb * (d - 1))) & (BB - 1)) - hB);
}

inline int64_t gcd(int64_t a, int64_t b) { a = abs(a); b = abs(b); while (b) { int64_t t = b; b = a % b; a = t; } return a; }

void trial_divide_side(mpz_t N, const vector<int>& fb_p, const vector<int>& small_primes, string& str, stringstream& stream) {
    if (small_primes.empty() || fb_p.empty()) return;
    int p = small_primes[0], k = 0, max_small = 1000;
    while (p < fb_p.back()) {
        while (mpz_fdiv_ui(N, p) == 0) {
            mpz_divexact_ui(N, N, p);
            if (!str.empty() && str.back() != ':' && str.back() != ',') str += ",";
            stream.str(""); stream << hex << p; 
            str += stream.str();
        }
        if (p < max_small) {
            if ((size_t)(k + 1) < small_primes.size()) p = small_primes[++k];
            else p = max_small + 1;
            if (p > max_small) { k = 0; while ((size_t)k < fb_p.size() && fb_p[k] < max_small) k++; }
        } else {
            if ((size_t)(k + 1) < fb_p.size()) p = fb_p[++k]; else break;
        }
    }
}

bool cofactorize_side(mpz_t N, mpz_t S, const mpz_class& lpb, string& str, int BASE,
                      stack<mpz_ptr>& QN, stack<int>& Q, int* algarr, mpz_t* pi,
                      mpz_t factor, mpz_t p1, mpz_t p2, mpz_t t)
{
    if (mpz_cmp_ui(N, 1) == 0) return true;
    if (mpz_probab_prime_p(N, 30)) {
        if (mpz_cmpabs(N, lpb.get_mpz_t()) > 0) return false;
        if (!str.empty() && str.back() != ':' && str.back() != ',') str += ",";
        char buf[32]; mpz_get_str(buf, BASE, N); str += buf; 
        return true;
    }
    while (!Q.empty()) Q.pop(); while (!QN.empty()) QN.pop();
    QN.push(N); Q.push(2); Q.push(1); Q.push(0); Q.push(3);
    while (!QN.empty()) {
        mpz_ptr Ncur = QN.top(); QN.pop();
        int l = Q.top(); Q.pop();
        int j = 0; bool factored = false;
        while (!factored) {
            int alg = Q.top(); Q.pop(); j++;
            if (alg == 0) factored = PollardPm1(Ncur, S, factor);
            else if (alg == 1) factored = EECM(Ncur, S, factor, 25921, 83521, 19, 9537, 2737);
            else if (alg == 2) factored = EECM(Ncur, S, factor, 1681, 707281, 3, 19642, 19803);
            if (!factored) { if (j >= l) return false; }
            else {
                mpz_set(p1, factor); mpz_divexact(p2, Ncur, factor);
                if (mpz_cmpabs(p1, p2) > 0) { mpz_set(t, p1); mpz_set(p1, p2); mpz_set(p2, t); }
                int lnext = l - j;
                for (int lt = 0; lt < lnext; lt++) algarr[lt] = Q.top(), Q.pop();
                for (int lt = 0; lt < lnext; lt++) Q.push(algarr[lnext - 1 - lt]);
                if (lnext) Q.push(lnext);
                
                if (mpz_probab_prime_p(p1, 30)) {
                    if (mpz_cmpabs(p1, lpb.get_mpz_t()) > 0) return false;
                    if (!str.empty() && str.back() != ':' && str.back() != ',') str += ",";
                    char buf[32]; mpz_get_str(buf, BASE, p1); str += buf;
                } else if (lnext) { mpz_set(pi[0], p1); QN.push(pi[0]); }
                else return false;
                
                if (mpz_probab_prime_p(p2, 30)) {
                    if (mpz_cmpabs(p2, lpb.get_mpz_t()) > 0) return false;
                    if (!str.empty() && str.back() != ':' && str.back() != ',') str += ",";
                    char buf[32]; mpz_get_str(buf, BASE, p2); str += buf;
                } else if (lnext) { mpz_set(pi[1], p2); QN.push(pi[1]); }
                else return false;
            }
        }
    }
    return true;
}

// ==================== CORE COFACTORIZATION ====================
void GetlcmScalar(int B, mpz_t S, int* primes, int nump) {
    mpz_t* tree = new mpz_t[nump];
    mpz_t pe, pe1;
    mpz_init(pe); mpz_init(pe1);

    int n = 0;
    int p = 2;
    while (p < B) {
        mpz_set_ui(pe, p);
        mpz_mul_ui(pe1, pe, p);
        while (mpz_cmp_ui(pe1, B) < 0) {
            mpz_set(pe, pe1);
            mpz_mul_ui(pe1, pe, p);
        }
        mpz_init(tree[n]);
        mpz_set(tree[n], pe);
        n++;
        p = primes[n];
    }
    mpz_clear(pe); mpz_clear(pe1);

    uint64_t treepos = n - 1;
    while (treepos > 0) {
        for (int i = 0; i <= treepos; i += 2) {
            if (i < treepos) mpz_lcm(tree[i/2], tree[i], tree[i + 1]);
            else mpz_set(tree[i/2], tree[i]);
        }
        for (int i = (treepos >> 1); i < treepos - 1; i++) mpz_set_ui(tree[i + 1], 1);
        treepos = treepos >> 1;
    }
    mpz_set(S, tree[0]);

    for (int i = 0; i < n; i++) mpz_clear(tree[i]);
    delete[] tree;
}

inline __int128 make_int128(uint64_t lo, uint64_t hi) {
    __int128 N = hi; N = N << 64; N += lo; return N;
}

inline __int128 gcd128(__int128 a, __int128 b) {
    a = a < 0 ? -a : a; b = b < 0 ? -b : b; __int128 c;
    while (b != 0) { c = b; b = a % c; a = c; }
    return a;
}

bool PollardPm1(mpz_t N, mpz_t S, mpz_t factor) {
    int bitlen = mpz_sizeinbase(N, 2);
    if (bitlen < 64) {
        __int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N,1));
        __int128 factor128 = 1;
        PollardPm1_int128(N128, S, factor128, 0, 0);
        mp_limb_t* fl = mpz_limbs_modify(factor, 2);
        fl[0] = factor128 & MASK64; fl[1] = factor128 >> 64;
        return (factor128 > 1 && factor128 < N128);
    }
    return PollardPm1_mpz(N, S, factor);
}

bool PollardPm1_mpz(mpz_t N, mpz_t S, mpz_t factor) {
    int L = mpz_sizeinbase(S, 2);
    mpz_t g; mpz_init_set_ui(g, 2);
    for (int i = 2; i <= L; i++) {
        mpz_mul(g, g, g); mpz_mod(g, g, N);
        if (mpz_tstbit(S, L - i) == 1) {
            mpz_mul_2exp(g, g, 1);
            if (mpz_cmpabs(g, N) >= 0) mpz_sub(g, g, N);
        }
    }
    mpz_sub_ui(g, g, 1); mpz_gcd(factor, N, g);
    bool result = mpz_cmpabs_ui(factor, 1) > 0 && mpz_cmpabs(factor, N) < 0;
    mpz_clear(g); return result;
}

bool PollardPm1_int128(__int128 N, mpz_t S, __int128 &factor, int64_t, int64_t) {
    int L = mpz_sizeinbase(S, 2); __int128 g = 2;
    for (int i = 2; i <= L; i++) {
        g = g * g % N;
        if (mpz_tstbit(S, L - i) == 1) {
            g = g * 2; if (g >= N) g -= N;
        }
    }
    g -= 1; factor = gcd128(N, g);
    return (factor > 1 && factor < N);
}

bool EECM(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0) {
    int bitlen = mpz_sizeinbase(N, 2);
    if (bitlen < 64) {
        __int128 N128 = make_int128(mpz_getlimbn(N,0), mpz_getlimbn(N,1));
        __int128 factor128 = 1;
        EECM_int128(N128, S, factor128, d, a, X0, Y0, Z0, 0, 0);
        mp_limb_t* fl = mpz_limbs_modify(factor, 2);
        fl[0] = factor128 & MASK64; fl[1] = factor128 >> 64;
        return (factor128 > 1 && factor128 < N128);
    }
    return EECM_mpz(N, S, factor, d, a, X0, Y0, Z0);
}

bool EECM_mpz(mpz_t N, mpz_t S, mpz_t factor, int d, int a, int X0, int Y0, int Z0) {
    mpz_t SX, SY, SZ, A, B, B2, B3, C, dC, B2mC, D, CaD, B2mCmD, E, EmD, F, AF, G, AG, aC, DmaC, H, Hx2, J;
    mpz_t X0aY0xB, X0aY0xB_mCmD, X, Y, Z, mulmod;
    mpz_init(SX); mpz_init_set_ui(SX, X0); mpz_init(SY); mpz_init_set_ui(SY, Y0); mpz_init(SZ); mpz_init_set_ui(SZ, Z0);
    mpz_init(A); mpz_init(B); mpz_init(B2); mpz_init(B3); mpz_init(C); mpz_init(dC); mpz_init(B2mC);
    mpz_init(D); mpz_init(CaD); mpz_init(B2mCmD); mpz_init(E); mpz_init(EmD); mpz_init(F); mpz_init(AF);
    mpz_init(G); mpz_init(AG); mpz_init(aC); mpz_init(DmaC); mpz_init(H); mpz_init(Hx2); mpz_init(J);
    mpz_init(X0aY0xB); mpz_init(X0aY0xB_mCmD); mpz_init(X); mpz_init(Y); mpz_init(Z); mpz_init(mulmod);

    int L = mpz_sizeinbase(S, 2);
    for(int i = 2; i <= L; i++) {
        mpz_add(B, SX, SY); mpz_mul(mulmod, B, B); mpz_mod(B2, mulmod, N);
        mpz_mul(mulmod, SX, SX); mpz_mod(C, mulmod, N);
        mpz_mul(mulmod, SY, SY); mpz_mod(D, mulmod, N);
        mpz_mul_ui(E, C, a); mpz_add(F, E, D);
        mpz_mul(mulmod, SZ, SZ); mpz_mod(H, mulmod, N);
        mpz_mul_2exp(Hx2, H, 1); mpz_sub(J, F, Hx2);
        mpz_add(CaD, C, D); mpz_sub(B2mCmD, B2, CaD);
        mpz_mul(X, B2mCmD, J); mpz_sub(EmD, E, D);
        mpz_mul(Y, F, EmD); mpz_mul(Z, F, J);
        mpz_mod(SX, X, N); mpz_mod(SY, Y, N); mpz_mod(SZ, Z, N);

        if(mpz_tstbit(S, L - i) == 1) {
            mpz_mul_ui(A, SZ, Z0); mpz_add(B, SX, SY);
            mpz_mul(mulmod, A, A); mpz_mod(B3, mulmod, N);
            mpz_mul_ui(C, SX, X0); mpz_mul_ui(D, SY, Y0);
            mpz_mul_ui(dC, C, d); mpz_add(CaD, C, D);
            mpz_mul(mulmod, dC, D); mpz_mod(E, mulmod, N);
            mpz_sub(F, B3, E); mpz_add(G, B3, E);
            mpz_mul_ui(mulmod, B, X0+Y0); mpz_mod(X0aY0xB, mulmod, N);
            mpz_sub(X0aY0xB_mCmD, X0aY0xB, CaD);
            mpz_mul(mulmod, A, F); mpz_mod(AF, mulmod, N);
            mpz_mul(X, AF, X0aY0xB_mCmD);
            mpz_mul(mulmod, A, G); mpz_mod(AG, mulmod, N);
            mpz_mul_ui(aC, C, a); mpz_sub(DmaC, D, aC);
            mpz_mul(Y, AG, DmaC); mpz_mul(Z, F, G);
            mpz_mod(SX, X, N); mpz_mod(SY, Y, N); mpz_mod(SZ, Z, N);
        }
    }
    mpz_gcd(factor, N, SX);

    mpz_clear(mulmod); mpz_clear(SZ); mpz_clear(SY); mpz_clear(SX);
    mpz_clear(A); mpz_clear(B); mpz_clear(B2); mpz_clear(B3); mpz_clear(C); mpz_clear(dC); mpz_clear(B2mC);
    mpz_clear(D); mpz_clear(CaD); mpz_clear(B2mCmD); mpz_clear(E); mpz_clear(EmD); mpz_clear(F); mpz_clear(AF);
    mpz_clear(G); mpz_clear(AG); mpz_clear(aC); mpz_clear(DmaC); mpz_clear(H); mpz_clear(Hx2); mpz_clear(J);
    mpz_clear(X0aY0xB); mpz_clear(X0aY0xB_mCmD); mpz_clear(X); mpz_clear(Y); mpz_clear(Z);

    return mpz_cmpabs_ui(factor, 1) > 0 && mpz_cmpabs(factor, N) < 0;
}

bool EECM_int128(__int128 N, mpz_t S, __int128 &factor, int d, int a, int X0, int Y0, int Z0, int64_t, int64_t) {
    __int128 SX = X0, SY = Y0, SZ = Z0;
    __int128 A, B, B2, B3, C, dC, B2mC, D, CaD, B2mCmD, E, EmD, F, AF, G, AG, aC, DmaC, H, Hx2, J;
    __int128 X0aY0xB, X0aY0xB_mCmD, X, Y, Z;

    int L = mpz_sizeinbase(S, 2);
    for(int i = 2; i <= L; i++) {
        B = SX + SY; B2 = (B * B) % N; C = (SX * SX) % N; D = (SY * SY) % N;
        E = C * a; F = E + D; H = (SZ * SZ) % N; Hx2 = H << 1; J = F - Hx2;
        CaD = C + D; B2mCmD = B2 - CaD; X = B2mCmD * J; EmD = E - D; Y = F * EmD; Z = F * J;
        SX = X % N; SY = Y % N; SZ = Z % N;

        if(mpz_tstbit(S, L - i) == 1) {
            A = SZ * Z0; B = SX + SY; B3 = (A * A) % N; C = SX * X0; D = SY * Y0;
            dC = C * d; CaD = C + D; E = (dC * D) % N; F = B3 - E; G = B3 + E;
            X0aY0xB = (B * (X0 + Y0)) % N; X0aY0xB_mCmD = X0aY0xB - CaD;
            AF = (A * F) % N; X = AF * X0aY0xB_mCmD; AG = (A * G) % N; aC = C * a; DmaC = D - aC;
            Y = AG * DmaC; Z = F * G; SX = X % N; SY = Y % N; SZ = Z % N;
        }
    }
    factor = gcd128(N, SX);
    return (factor > 1 && factor < N);
}

// ==================== MAIN ====================
int main(int argc, char** argv)
{
    MASK64 = ((int128_t)1 << 64) - 1;
    if (argc != 17) {
        cerr << "Usage: ./hdlcs inputpoly fb d Amax Bmax N Bmin Bmax Rmin Rmax th0 th1 lpb bits mbb bb" << endl;
        return 1;
    }

    int d = atoi(argv[3]);
    mpz_class maxA(argv[4]), maxB(argv[5]), lpb(argv[13]);
    int N_units = atoi(argv[6]), Bmin = atoi(argv[7]), Bmax = atoi(argv[8]);
    int Rmin = atoi(argv[9]), Rmax = atoi(argv[10]);
    uint8_t th[2] = {(uint8_t)atoi(argv[11]), (uint8_t)atoi(argv[12])};
    int64_t cofmax = 1LL << atoi(argv[14]);
    int mbb = atoi(argv[15]), bb = atoi(argv[16]);

    mpz_poly f0, f1;
    mpz_poly_init(f0, 10); mpz_poly_init(f1, 10);
    double skew = 1.0; int degf, degg;
    parse_polynomial(argv[1], f0, f1, skew, degf, degg);

    SieveSide sides[2];
    load_factor_base(argv[2], sides, th);

    vector<int> small_primes;
    {
        int maxs = 1 << 21; vector<char> sieve(maxs + 1, 0);
        for (int i = 2; i * i <= maxs; i++) if (!sieve[i]) for (int j = i * i; j <= maxs; j += i) sieve[j] = 1;
        for (int i = 2; i <= maxs; i++) if (!sieve[i]) small_primes.push_back(i);
    }

    gmp_randstate_t state; gmp_randinit_default(state); gmp_randseed_ui(state, 123ul);
    vector<mpz_class> Ai(d - 2), Bi(d - 2);
    vector<mpz_ptr> Ai_ptr(d - 2), Bi_ptr(d - 2);

    for (int i = 0; i < d - 2; i++) { 
        Ai_ptr[i] = Ai[i].get_mpz_t(); 
        Bi_ptr[i] = Bi[i].get_mpz_t(); 
    }

    mpz_poly i1; mpz_poly_init(i1, 3);
    mpz_t N0, N1, S, factor, p1, p2, t, g1, A, B, pi[8];
    mpz_init(N0); mpz_init(N1); mpz_init(S); mpz_init(factor); mpz_init(p1); mpz_init(p2);
    mpz_init(t); mpz_init(g1); mpz_init(A); mpz_init(B);
    for (int i = 0; i < 8; i++) mpz_init(pi[i]);
    GetlcmScalar(cofmax, S, small_primes.data(), small_primes.size());

    // Allocate on the heap
    std::unique_ptr<MITM_Workspace> alg_ws_ptr(new MITM_Workspace());
    std::unique_ptr<MITM_Workspace> rat_ws_ptr(new MITM_Workspace());

    // Create references so the rest of your code (thread_func, etc.) doesn't need to change
    MITM_Workspace& alg_ws = *alg_ws_ptr;
    MITM_Workspace& rat_ws = *rat_ws_ptr;

    int8_t alpha[4] = {-2, -1, 0, 1};

    // --- Thread Function designed to be pinned to Core/Hyperthread ---
    auto thread_func = [&](int side_idx, int thread_id) {
        MITM_Workspace& ws = (side_idx == 0) ? alg_ws : rat_ws;
        const SieveSide& side = sides[side_idx];
        int kmax = side.k;
        uint8_t last_logp = 0;
        
        for (int i = 0; i < kmax; i++) {
            if (side.p[i] < Bmin) continue;
            if (side.p[i] >= Bmax) break;
            
            int64_t p = side.p[i];
            uint8_t logp = (uint8_t)max(1.0, log2(p));
            
            // OPTIMIZATION: Trigger global thread flush if logp has advanced
            if (logp != last_logp) {
                if (last_logp != 0) {
                    ws.flush_all_buckets(thread_id, last_logp);
                }
                last_logp = logp;
            }

            int ni = side.n[i];
            
            for (int j = 0; j < ni; j++) {
                int r = side.r[side.r_offset[i] + j];
                
                // Map logical matrix variables (rel2A indexing) into physical loop variables
                int64_t a[16];
                a[0] = p; // Dim 15 (mapped from bits 30..31)
                for (int k = 0; k < 14; k++) {
                    a[k + 1] = (mpz_fdiv_ui(Ai_ptr[k], p) * r + mpz_fdiv_ui(Bi_ptr[k], p)) % p;
                } // Dims 0..13 (mapped from bits 0..27)
                a[15] = r; // Dim 14 (mapped from bits 28..29)
                
                int64_t v[16][4];
                for (int m = 0; m < 16; m++) {
                    for (int k = 0; k < 4; k++) v[m][k] = (int64_t)alpha[k] * a[m];
                }
                
                ws.current_run++;
                
                // Phase 1: Thread 0 populates the hash table with Left Half
                if (thread_id == 0) {
                    ws.node_count = 1;
                    int64_t s1 = 0;
                    for (int m = 1; m <= 7; m++) s1 += v[m][0];
                    
                    for (uint32_t left_idx = 0; left_idx < 16384; left_idx++) {
                        int64_t m1 = s1 % p; if (m1 < 0) m1 += p;
                        uint32_t h = hash_func(m1);
                        
                        while (ws.last_run[h] == ws.current_run && ws.keys[h] != m1) h = (h + 1) & HASH_MASK;
                        if (ws.last_run[h] != ws.current_run) {
                            ws.last_run[h] = ws.current_run;
                            ws.keys[h] = m1;
                            ws.head[h] = 0;
                        }
                        
                        uint32_t node = ws.node_count++;
                        ws.nodes[node].coords = left_idx;
                        ws.nodes[node].next = ws.head[h];
                        ws.head[h] = node;
                        
                        // Advance Left Odometer
                        for (int m = 0; m < 7; m++) {
                            int idx_v = m + 1;
                            int curr_val = (left_idx >> (2 * m)) & 3;
                            if (curr_val < 3) {
                                s1 += (v[idx_v][curr_val + 1] - v[idx_v][curr_val]);
                                break;
                            } else {
                                s1 += (v[idx_v][0] - v[idx_v][3]);
                            }
                        }
                    }
                }
                
                // Block until Thread 0 finishes caching Left Half
                ws.barrier.wait();
                
                // Phase 2: Parallel execution of Right Half
#ifdef DEBUG
                // For single thread debug, thread_id 0 covers the full range
                uint32_t start_idx = 0;
                uint32_t end_idx = 65536;
#else
                uint32_t start_idx = (thread_id == 0) ? 0 : 32768;
                uint32_t end_idx = (thread_id == 0) ? 32768 : 65536;
#endif
                
                int64_t s2 = 0;
                for (int m = 0; m < 8; m++) {
                    s2 += v[m + 8][(start_idx >> (2 * m)) & 3];
                }
                
                for (uint32_t right_idx = start_idx; right_idx < end_idx; right_idx++) {
                    int64_t m2 = s2 % p; if (m2 < 0) m2 += p;
                    int64_t target = p - m2; if (target == p) target = 0;
                    
                    uint32_t h = hash_func(target);
                    while (ws.last_run[h] == ws.current_run) {
                        if (ws.keys[h] == target) {
                            uint32_t curr = ws.head[h];
                            while (curr != 0) {
                                uint32_t left_idx = ws.nodes[curr].coords;
                                // Shift variables into correct mapping for rel2A logic
                                uint32_t base_idx = left_idx | (right_idx << 14);
                                for (int k = 0; k < 4; k++) {
                                    uint32_t full_idx = base_idx | (k << 30);
                                    ws.emit_to_bucket(thread_id, full_idx, logp);
                                }
                                curr = ws.nodes[curr].next;
                            }
                            break;
                        }
                        h = (h + 1) & HASH_MASK;
                    }
                    
                    // Advance Right Odometer
                    for (int m = 0; m < 8; m++) {
                        int idx_v = m + 8;
                        int curr_val = (right_idx >> (2 * m)) & 3;
                        if (curr_val < 3) {
                            s2 += (v[idx_v][curr_val + 1] - v[idx_v][curr_val]);
                            break;
                        } else {
                            s2 += (v[idx_v][0] - v[idx_v][3]);
                        }
                    }
                }
                
                // Synchronize before tearing down or starting next prime
                ws.barrier.wait();
            }
        }

        // Final cleanup for any residual buffer values
        if (last_logp != 0) {
            ws.flush_all_buckets(thread_id, last_logp);
        }
    };

    for (int nn = 0; nn < N_units; nn++) {
        cout << "\n# ========== Unit " << (nn + 1) << " / " << N_units << " ==========" << endl;
        
        for (int i = 0; i < d - 2; i++) {
            mpz_urandomm(Ai_ptr[i], state, maxA.get_mpz_t());
            mpz_urandomm(Bi_ptr[i], state, maxB.get_mpz_t());
            char str1[64], str2[64];
            mpz_get_str(str1, 10, Ai_ptr[i]); mpz_get_str(str2, 10, Bi_ptr[i]);
            cout << "# " << str1 << "*x + " << str2 << endl;
        }

        // Clear 4GB sieve arrays
        memset(alg_ws.array, 0, ARRAY_SIZE);
        memset(rat_ws.array, 0, ARRAY_SIZE);

        auto start = std::chrono::high_resolution_clock::now();
        
#ifdef DEBUG
        cout << "# DEBUG: Launching serial execution for HDLCS MITM..." << endl;
        thread_func(0, 0); // Alg Side
        thread_func(1, 0); // Rat Side
#else
        cout << "# Launching 4-thread execution for HDLCS parallel MITM..." << endl;
        // Exactly 4 threads (2 Alg, 2 Rat)
        std::thread t0(thread_func, 0, 0);
        std::thread t1(thread_func, 0, 1);
        std::thread t2(thread_func, 1, 0);
        std::thread t3(thread_func, 1, 1);
        
        t0.join(); t1.join(); t2.join(); t3.join();
#endif
        
        auto end = std::chrono::high_resolution_clock::now();
        cout << "# HDLCS Enumeration Finished! Time: " 
             << std::chrono::duration<double>(end - start).count() << "s" << endl;

        // --- Linear Sweep utilizing the same 4 threads ---
        start = std::chrono::high_resolution_clock::now();
        
        vector<int64_t> common;
        std::mutex common_mutex;
        
        auto sweep_worker = [&](uint64_t start_idx, uint64_t end_idx) {
            vector<int64_t> local_common;
            for (uint64_t idx = start_idx; idx < end_idx; idx++) {
                if (alg_ws.array[idx] >= th[0] && rat_ws.array[idx] >= th[1]) {
                    int64_t A64 = rel2A(d, Ai_ptr.data(), idx, bb);
                    int64_t B64 = rel2B(d, Bi_ptr.data(), idx, bb);
                    int64_t g_val = gcd(A64, B64);
                    if (A64 != 0 && B64 != 0 && abs(A64 / g_val) != 1) {
                        local_common.push_back(idx);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(common_mutex);
            common.insert(common.end(), local_common.begin(), local_common.end());
        };

#ifdef DEBUG
        cout << "# DEBUG: Executing serial linear sweep..." << endl;
        sweep_worker(0, ARRAY_SIZE);
#else
        cout << "# Executing concurrent multi-threaded linear sweep..." << endl;
        uint64_t chunk_size = ARRAY_SIZE / 4;
        std::thread s0(sweep_worker, 0, chunk_size);
        std::thread s1(sweep_worker, chunk_size, 2 * chunk_size);
        std::thread s2(sweep_worker, 2 * chunk_size, 3 * chunk_size);
        std::thread s3(sweep_worker, 3 * chunk_size, ARRAY_SIZE);
        
        s0.join(); s1.join(); s2.join(); s3.join();
#endif
        
        end = std::chrono::high_resolution_clock::now();
        cout << "# Linear Sweep Finished! Time: " 
             << std::chrono::duration<double>(end - start).count() << "s" << endl;
        cout << "# " << common.size() << " potential relations found." << endl;

        // --- Cofactorization ---
        cout << "# Starting cofactorization..." << endl;
        start = std::chrono::high_resolution_clock::now();
        int R = 0;
        int samples = 0;
        stringstream stream;
        
        for (int64_t id : common) {
            mpz_set_si(A, rel2A(d, Ai_ptr.data(), id, bb));
            mpz_set_si(B, rel2B(d, Bi_ptr.data(), id, bb));
            mpz_gcd(g1, A, B); mpz_divexact(A, A, g1); mpz_divexact(B, B, g1);
            mpz_poly_setcoeff(i1, 1, A); mpz_poly_setcoeff(i1, 0, B);
            mpz_poly_resultant(N0, f0, i1); mpz_poly_resultant(N1, f1, i1);
            mpz_abs(N0, N0); mpz_abs(N1, N1);

            string str = mpz_get_str(NULL, 10, A) + (string)"," + mpz_get_str(NULL, 10, B) + ":";
            stack<mpz_ptr> QN; stack<int> Q; int algarr[3];

            trial_divide_side(N0, sides[0].p, small_primes, str, stream);
            if (cofactorize_side(N0, S, lpb, str, 16, QN, Q, algarr, pi, factor, p1, p2, t)) {
                str += ":";
                trial_divide_side(N1, sides[1].p, small_primes, str, stream);
                if (cofactorize_side(N1, S, lpb, str, 16, QN, Q, algarr, pi, factor, p1, p2, t)) {
                    cout << str << endl;
                    R++;
                }
            }
            if (samples < 10) {
                cout << str << endl;
                samples++;
            }
        }
        end = std::chrono::high_resolution_clock::now();
        cout << "# Cofactorization took " << std::chrono::duration<double>(end - start).count() << "s" << endl;
        cout << "# " << R << " actual relations found." << endl;
    }

    // Cleanup
    mpz_poly_clear(f0); mpz_poly_clear(f1); mpz_poly_clear(i1);
    mpz_clear(N0); mpz_clear(N1); mpz_clear(S); mpz_clear(factor);
    for (int i = 0; i < 8; i++) mpz_clear(pi[i]);
    gmp_randclear(state);
    return 0;
}

