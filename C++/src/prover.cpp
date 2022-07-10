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
    for(int i = 0; i < irow; i++) {
        result[i] = inner_productp(vec, input[i], vrow);
    }
    return result;
}

Proof fliop(uint64_t** input, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    uint64_t L = var;
    uint64_t T = copy;

    uint64_t s = (T - 1) / k + 1;

    uint64_t eta = generate_challenge();

    //Prepare Input
    begin_time = clock();

    uint64_t eta_power = 1;
    uint64_t* meta_left = new uint64_t[2 * s * k];
    uint64_t** input_left = new uint64_t*[2 * s];
    for(int i = 0; i < s; i++) {
        input_left[2*i] = meta_left + 2 * i * k;
        input_left[2*i+1] = meta_left + 2 * i * k + k;
        for(int j = 0; j < k; j++) {
            if(i * k + j >= T) {
                input_left[2*i][j] = 0;
                input_left[2*i+1][j] = 0;
            }
            else {
                input_left[2*i][j] = mul_modp(input[0][i * k + j], eta_power);
                input_left[2*i+1][j] = mul_modp(input[2][i * k + j], eta_power);
                eta_power = mul_modp(eta_power, eta);
            }
        }
    }
    
    uint64_t* meta_right = new uint64_t[2 * s * k];
    uint64_t** input_right = new uint64_t*[2 * s];
    for(int i = 0; i < s; i++) {
        input_right[2*i] = meta_right + 2 * i * k;
        input_right[2*i+1] = meta_right + 2 * i * k + k;
        for(int j = 0; j < k; j++) {
            if(i * k + j >= T) {
                input_right[2 * i][j] = 0;
                input_right[2 * i + 1][j] = 0;
            }
            else {
                input_right[2 * i][j] = input[1][i * k + j];
                input_right[2 * i + 1][j] = input[3][i * k + j];
            }
        }
    }
    finish_time = clock();
    cout<<"Prepare Input Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    s *= 2;

    vector< vector<uint64_t> > p_coeffs_ss1;
    vector< vector<uint64_t> > p_coeffs_ss2;

    uint64_t** base = get_bases(k);

    while(true){
        cout<<"s : "<<s<<endl;
        cout<<"k : "<<k<<endl;

        //interpolation
        begin_time = clock();
        uint64_t** eval_left_polys = new uint64_t*[2*k-1];
        for(int i = 0; i < k; i++){
            eval_left_polys[i] = new uint64_t[s];
            for(int j = 0; j < s; j++){
                eval_left_polys[i][j] = input_left[j][i];
            }
        }
        for(int i = k; i < 2 * k - 1; i++){
            eval_left_polys[i] = matrix_mul(base[i - k], input_left, k, s);
        }
        uint64_t** eval_right_polys = new uint64_t*[2*k-1];
        for(int i = 0; i < 2 * k - 1; i++){
            if(i < k) {
                eval_right_polys[i] = new uint64_t[s];
                for(int j = 0; j < s; j++){
                    eval_right_polys[i][j] = input_right[j][i];
                }
            }
            else {
                eval_right_polys[i] = matrix_mul(base[i - k], input_right, k, s);
            }
        }
        finish_time = clock();
        cout<<"Interpolation Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

        //compute p(X)
        begin_time = clock();
        uint64_t* eval_p_poly = new uint64_t[2 * k - 1];
        for(int i = 0; i < 2 * k - 1; i++) {
            eval_p_poly[i] = inner_productp(eval_left_polys[i], eval_right_polys[i], s);
        }
        finish_time = clock();
        cout<<"Compute P(X) Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

        //generate proof
        begin_time = clock();
        vector<uint64_t> ss1, ss2;
        for(int i = 0; i < 2 * k - 1; i++) {
            ss1.push_back(get_rand());
            ss2.push_back(eval_p_poly[i] - ss1[i]);
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
        uint64_t* eval_base = evaluate_bases(k, r);
        for(int i = 0; i < s; i++) {
            input_left[i][0] = inner_productp(eval_base, input_left[i], k);
            input_right[i][0] = inner_productp(eval_base, input_right[i], k);
        }
        uint64_t s0 = s;
        s = (s - 1) / k + 1;
        for(int i = 0; i < s; i++) {
            for(int j = 0; j < k; j++) {
                if (i * k + j < s0) {
                    input_left[i][j] = input_left[i * k + j][0];
                    input_right[i][j] = input_right[i * k + j][0];
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

Proof prove_and_gate(uint64_t _party_id, uint64_t** inputs, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid) {
    return fliop(inputs, var, copy, k, sid);
}

int main() {
    uint64_t T = 10000000;
    uint64_t L = 5;
    uint64_t k = 10;
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

