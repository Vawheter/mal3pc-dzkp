#include<cstdio>
#include<cmath>
#include<iostream>

typedef unsigned __int128 uint128_t;

using namespace std;

static const uint32_t PRIME_EXP = 61;
static const uint64_t PR = 2305843009213693951;

uint64_t modp(uint64_t a) {
    uint64_t res = (a>>PRIME_EXP) + (a & PR);
    if (res >= PR) {
        res -= PR;
    }
    return res;
}

uint64_t modp_128(uint128_t a){
    uint64_t higher, middle, lower;
    higher = (a >> (2 * PRIME_EXP));
    middle = (a >> PRIME_EXP) & PR;
    lower = a & PR;
    return modp(higher + middle + lower);
}

uint64_t neg_modp(uint64_t a) {
    return PR - a;
}

uint64_t add_modp(uint64_t a, uint64_t b) {
    uint64_t res = a + b;
    if (res >= PR) {
        res -= PR;
    }
    return res;
}

uint64_t mul_modp(uint64_t a, uint64_t b) {
    uint128_t res = ((uint128_t) a) * ((uint128_t) b);
    uint64_t higher = (res>>PRIME_EXP);
    uint64_t lower = res & PR;
    return add_modp(higher, lower);
}

uint64_t inverse(uint64_t a) {
    uint64_t left = a;
    uint64_t right = PR;
    uint64_t x = 1, y = 0, u = 0, v = 1;
    uint64_t gcd = a;
    uint64_t w, z;
    while(left != 0) {
        w = right / left;
        z = right % left;
        right = left;
        left = z;

        z = u - w * x;
        u = x;
        x = z;

        z = v - w * y;
        v = y;
        y = z;
    }
    if (u >= PR) {
        u += PR;
    }
    return u;
}

uint64_t inner_productp(uint64_t* a, uint64_t* b, uint64_t size) {
    uint128_t result = 0;
    uint64_t bound = 63;
    uint64_t start, end;
    start = 0;
    while(true) {
        if (start + bound < size) {
            end = start + bound;
        }
        else {
            end = size;
        }
        for(int i = start; i < end; i++) {
            result += ((uint128_t)a[i]) * ((uint128_t)b[i]);
        }
        result = modp_128(result);
        start = end;
        if (start == size) break;
    }
    return result;
}