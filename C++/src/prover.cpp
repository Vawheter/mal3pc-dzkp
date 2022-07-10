#include "arithmetic.h"
#include <cstdlib>
#include <ctime>
#include <vector>

clock_t begin_time, finish_time;

uint64_t get_rand() {
    srand((unsigned)time(NULL));
    return rand() & PR;
}

uint64_t generate_challenge() {
    return get_rand();
}

struct Proof {
    vector< vector<uint64_t> > p_coeffs_ss1;
    vector< vector<uint64_t> > p_coeffs_ss2;
};

uint64_t** get_bases(uint64_t n) {
    uint64_t** result = new uint64_t*[n-1];
    for (int i = 0; i < n - 1; i++) {
        result[i] = new uint64_t[n];
        for(int j = 0; j < n; j++) {
            result[i][j] = 1;
            for(int l = 0; l < n; l++) {
                if (l != j) {
                    uint64_t denominator, numerator;
                    if (j > l) {
                        denominator = j - l;
                    }
                    else {
                        denominator = neg_modp(l - j);
                    }
                    numerator = i + n - l;
                    result[i][j] = mul_modp(result[i][j], mul_modp(inverse(denominator), numerator));
                }
            }
        }
    }
    return result;
}

uint64_t* evaluate_bases(uint64_t n, uint64_t r) {
    uint64_t* result = new uint64_t[n];
    for(int i = 0; i < n; i++) {
        result[i] = 1;
        for(int j = 0; j < n; j++) {
            if (j != i) {
                uint64_t denominator, numerator; 
                if (i > j) { 
                    denominator = i - j;
                } 
                else { 
                    denominator = neg_modp(j - i);
                }
                if (r > j) { 
                    numerator = r - j; 
                } 
                else { 
                    numerator = neg_modp(j - r);
                }
                result[i] = mul_modp(result[i], mul_modp(inverse(denominator), numerator));
            }
        }
    }
    return result;
}

uint64_t* matrix_mul(uint64_t* vec, uint64_t** input, uint64_t vrow, uint64_t irow) {
    uint64_t* result = new uint64_t[irow];
    uint128_t temp_result;
    uint64_t higher, middle, lower;
    for(int i = 0; i < irow; i++) {
        temp_result = 0;
        for(int j = 0; j < vrow; j++) {
            temp_result += ((uint128_t) vec[j]) * ((uint128_t) input[j][i]);
        }
        higher = (temp_result >> (2 * PRIME_EXP));
        middle = (temp_result >> PRIME_EXP) & PR;
        lower = temp_result & PR;
        result[i] = modp(higher + middle + lower);
    }

    return result;
}

Proof fliop(uint64_t** input, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    uint64_t L = var;
    uint64_t T = copy;

    uint64_t s = (T - 1) / k + 1;

    uint64_t eta = generate_challenge();

    //Prepare Another Input
    begin_time = clock();

    uint64_t eta_power = 1;
    uint64_t* meta_left_prime = new uint64_t[2 * s * k];
    uint64_t** input_left_prime = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        input_left_prime[i] = meta_left_prime + i * 2 * s;
        for(int j = 0; j < s; j++) {
            if(i * s + j >= T) {
                input_left_prime[i][2 * j] = 0;
                input_left_prime[i][2 * j + 1] = 0;
            }
            else {
                input_left_prime[i][2 * j] = mul_modp(input[0][i * s + j], eta_power);
                input_left_prime[i][2 * j + 1] = mul_modp(input[2][i * s + j], eta_power);
                eta_power = mul_modp(eta_power, eta);
            }
        }
    }

    uint64_t* meta_right_prime = new uint64_t[2 * s * k];
    uint64_t** input_right_prime = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        input_right_prime[i] = meta_right_prime + i * 2 * s;
        for(int j = 0; j < s; j++) {
            if(i * s + j >= T) {
                input_right_prime[i][2 * j] = 0;
                input_right_prime[i][2 * j + 1] = 0;
            }
            else {
                input_right_prime[i][2 * j] = input[1][i * s + j];
                input_right_prime[i][2 * j + 1] = input[3][i * s + j];
            }
        }
    }
    finish_time = clock();
    cout<<"Prepare Input Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    s *= 2;

    vector< vector<uint64_t> > p_coeffs_ss1;
    vector< vector<uint64_t> > p_coeffs_ss2;

    uint64_t** base = get_bases(k);
    uint64_t* eval_base;
    uint64_t s0;

    while(true){
        cout<<"s : "<<s<<endl;
        cout<<"k : "<<k<<endl;

        //interpolation
        begin_time = clock();
        uint64_t** eval_left_polys = new uint64_t*[k-1];
        for(int i = 0; i < k - 1; i++){
            eval_left_polys[i] = matrix_mul(base[i], input_left_prime, k, s);
        }
        uint64_t** eval_right_polys = new uint64_t*[k-1];
        for(int i = 0; i < k - 1; i++){
            eval_right_polys[i] = matrix_mul(base[i], input_right_prime, k, s);
        }
        finish_time = clock();
        cout<<"Interpolation Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

        //compute p(X)
        begin_time = clock();
        uint64_t* eval_p_poly = new uint64_t[2 * k - 1];
        for(int i = 0; i < k; i++) {
            eval_p_poly[i] = inner_productp(input_left_prime[i], input_right_prime[i], s);
        }
        for(int i = 0; i < k - 1; i++) {
            eval_p_poly[i + k] = inner_productp(eval_left_polys[i], eval_right_polys[i], s);
        }
        finish_time = clock();
        cout<<"Compute P(X) Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

        //generate proof
        begin_time = clock();
        vector<uint64_t> ss1, ss2;
        uint64_t temp;
        for(int i = 0; i < 2 * k - 1; i++) {
            ss1.push_back(get_rand());
            if(eval_p_poly[i] > ss1[i]) {
                temp = eval_p_poly[i] - ss1[i];
            }
            else {
                temp = PR - ss1[i] + eval_p_poly[i];
            }
            ss2.push_back(temp);
        }
        p_coeffs_ss1.push_back(ss1);
        p_coeffs_ss2.push_back(ss2);
        finish_time = clock();
        cout<<"Generate Proof Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

        if (s == 1) {
            break;
        }

        // Prepare Next Input
        begin_time = clock();
        uint64_t r = generate_challenge();
        eval_base = evaluate_bases(k, r);
        s0 = s;
        s = (s - 1) / k + 1;
        uint128_t temp_result;
        uint64_t higher, middle, lower, index;
        for(int i = 0; i < k; i++) {
            for(int j = 0; j < s; j++) {
                index = i * s + j;
                if (index < s0) {
                    temp_result = 0;
                    for(int l = 0; l < k; l++) {
                        temp_result += ((uint128_t) eval_base[l]) * ((uint128_t) input_left_prime[l][index]);
                    }
                    higher = (temp_result >> (2 * PRIME_EXP));
                    middle = (temp_result >> PRIME_EXP) & PR;
                    lower = temp_result & PR;
                    input_left_prime[i][j] = modp(higher + middle + lower);

                    temp_result = 0;
                    for(int l = 0; l < k; l++) {
                        temp_result += ((uint128_t) eval_base[l]) * ((uint128_t) input_right_prime[l][index]);
                    }
                    higher = (temp_result >> (2 * PRIME_EXP));
                    middle = (temp_result >> PRIME_EXP) & PR;
                    lower = temp_result & PR;
                    input_right_prime[i][j] = modp(higher + middle + lower);
                }
                else {
                    input_left_prime[i][j] = 0;
                    input_right_prime[i][j] = 0;
                }
            }
        }
        finish_time = clock();
        cout<<"Prepare Input Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;
    }

    Proof result = {p_coeffs_ss1, p_coeffs_ss2};
    return result;
}

Proof prove_and_gate(uint64_t _party_id, uint64_t** inputs, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    return fliop(inputs, var, copy, k, sid);
}

int main() {
    uint64_t T = 1000000;
    uint64_t L = 5;
    uint64_t k = 4;
    uint64_t _party_id = 1;

    cout<<"T: "<<T<<endl;
    uint64_t sid = get_rand();

    uint64_t** input = new uint64_t*[L];
    for(int i = 0; i < L; i++) {
        input[i] = new uint64_t[T];
        for(int j = 0; j < T; j++) {
            input[i][j] = get_rand();
        }
    }

    clock_t start, end;

    start = clock();
    Proof proof = prove_and_gate(_party_id, input, L, T, k, sid);
    end = clock();
    cout<<"Total Time = "<<double(end-start)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    return 0;
}

