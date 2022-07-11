#include "arithmetic.h"
#include <cstdlib>
#include <ctime>
#include <vector>

clock_t begin_time, finish_time;

uint64_t get_rand() {
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

Proof fliop(uint64_t** input_left, uint64_t** input_right, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    uint64_t L = var;
    uint64_t T = copy;
    uint64_t s = (T - 1) / k + 1;
    uint64_t eta = generate_challenge();

    //Prepare Input
    begin_time = clock();
    uint64_t eta_power = 1;
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < s; j++) {
            input_left[i][2 * j] = mul_modp(input_left[i][2 * j], eta_power);
            input_left[i][2 * j + 1] = mul_modp(input_left[i][2 * j + 1], eta_power);
            eta_power = mul_modp(eta_power, eta);
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
    uint64_t** eval_result = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        eval_result[i] = new uint64_t[k];
    }
    uint64_t* eval_p_poly = new uint64_t[2 * k - 1];
    uint64_t r;
    uint128_t temp_result;
    uint64_t index;

    while(true){
        cout<<"s : "<<s<<endl;
        cout<<"k : "<<k<<endl;

        //Compute P(X)
        begin_time = clock();
        for(int i = 0; i < k; i++) {
            for(int j = 0; j < k; j++) {
                eval_result[i][j] = inner_productp(input_left[i], input_right[j], s);
            }
        }

        for(int i = 0; i < k; i++) {
            eval_p_poly[i] = eval_result[i][i];
        }
        for(int i = 0; i < k - 1; i++) {
            eval_p_poly[i + k] = 0;
            for(int j = 0; j < k; j++) {
                for (int l = 0; l < k; l++) {
                    eval_p_poly[i+k] = add_modp(eval_p_poly[i+k], mul_modp(base[i][j], mul_modp(eval_result[j][l], base[i][l])));
                }
            }
        }
        finish_time = clock();
        cout<<"Interpolation Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

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
        r = generate_challenge();
        eval_base = evaluate_bases(k, r);
        s0 = s;
        s = (s - 1) / k + 1;
        for(int i = 0; i < k; i++) {
            for(int j = 0; j < s; j++) {
                index = i * s + j;
                if (index < s0) {
                    temp_result = 0;
                    for(int l = 0; l < k; l++) {
                        temp_result += ((uint128_t) eval_base[l]) * ((uint128_t) input_left[l][index]);
                    }
                    input_left[i][j] = modp_128(temp_result);

                    temp_result = 0;
                    for(int l = 0; l < k; l++) {
                        temp_result += ((uint128_t) eval_base[l]) * ((uint128_t) input_right[l][index]);
                    }
                    input_right[i][j] = modp_128(temp_result);
                }
                else {
                    input_left[i][j] = 0;
                    input_right[i][j] = 0;
                }
            }
        }
        finish_time = clock();
        cout<<"Prepare Input Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;
    }

    Proof result = {p_coeffs_ss1, p_coeffs_ss2};
    return result;
}

Proof prove_and_gate(uint64_t _party_id, uint64_t** input_left, uint64_t** input_right, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    return fliop(input_left, input_right, var, copy, k, sid);
}

void shape(uint64_t** input, uint64_t L, uint64_t T, uint64_t k, uint64_t** &input_left, uint64_t** &input_right) {
    uint64_t s = (T - 1) / k + 1;
    uint64_t* meta_left = new uint64_t[2 * s * k];
    input_left = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        input_left[i] = meta_left + i * 2 * s;
        for(int j = 0; j < s; j++) {
            if(i * s + j >= T) {
                input_left[i][2 * j] = 0;
                input_left[i][2 * j + 1] = 0;
            }
            else {
                input_left[i][2 * j] = input[0][i * s + j];
                input_left[i][2 * j + 1] = input[2][i * s + j];
            }
        }
    }

    uint64_t* meta_right = new uint64_t[2 * s * k];
    input_right = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        input_right[i] = meta_right + i * 2 * s;
        for(int j = 0; j < s; j++) {
            if(i * s + j >= T) {
                input_right[i][2 * j] = 0;
                input_right[i][2 * j + 1] = 0;
            }
            else {
                input_right[i][2 * j] = input[1][i * s + j];
                input_right[i][2 * j + 1] = input[3][i * s + j];
            }
        }
    }
}

int main() {
    uint64_t T = 10000000;
    uint64_t L = 5;
    uint64_t k = 4;
    uint64_t _party_id = 1;
    srand((unsigned)time(NULL));
    uint64_t sid = get_rand();

    uint64_t** input = new uint64_t*[L];
    for(int i = 0; i < L; i++) {
        input[i] = new uint64_t[T];
        for(int j = 0; j < T; j++) {
            input[i][j] = get_rand();
        }
    }

    uint64_t** input_left, **input_right;

    shape(input, L, T, k, input_left, input_right);

    cout<<"T: "<<T<<endl;
    clock_t start, end;

    start = clock();
    Proof proof = prove_and_gate(_party_id, input_left, input_right, L, T, k, sid);
    end = clock();
    cout<<"Total Time = "<<double(end-start)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    return 0;
}

