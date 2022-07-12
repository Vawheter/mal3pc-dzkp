#include "arithmetic.h"
#include <cstdlib>
#include <ctime>
#include <vector>

using namespace std;

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

bool test_evaluate_bases() {
    uint64_t k = 10;
    for(int i = 0; i < k; i++) {
        uint64_t* eval_base = evaluate_bases(k, i);
        if(eval_base[i] != 1) {
            cout << "evaluate_bases() incorrect" << endl;
            return false;
        }
        for(int j = 0; j < k; j++) {
            if(j != i) {
                if(eval_base[j] != 0) {
                    cout << "evaluate_bases() incorrect" << endl;
                    return false;
                }
            }
        }
    }
    cout << "evaluate_bases() correct" << endl;
    return true;
}


Proof fliop(uint64_t** input_left, uint64_t** input_right, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid, uint64_t* rands) {
    uint64_t L = var;
    uint64_t T = copy;
    uint64_t s = (T - 1) / k + 1;
    // uint64_t eta = generate_challenge();
    uint64_t eta = rands[0];

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

    uint16_t cnt = 0;

    while(true){
        cnt++;
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

        // Check p's evaluations
        for(int l = 0; l < 2 * k - 1; l++) {
            uint64_t* eval_base = evaluate_bases(k, l);
            uint64_t a = 0, b = 0, c = 0, d = 0, res = 0;
            for(int j = 0; j < s; j++) {
                uint128_t tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0;
                for(int i = 0; i < k; i++) {
                    tmp1 += ((uint128_t) eval_base[i]) * ((uint128_t) input_left[i][2 * j]);
                    tmp2 += ((uint128_t) eval_base[i]) * ((uint128_t) input_right[i][2 * j]);
                    tmp3 += ((uint128_t) eval_base[i]) * ((uint128_t) input_left[i][2 * j + 1]);
                    tmp4 += ((uint128_t) eval_base[i]) * ((uint128_t) input_right[i][2 * j + 1]);
                }
                a = modp_128(tmp1);
                b = modp_128(tmp2);
                c = modp_128(tmp3);
                d = modp_128(tmp4);
                if(l < k) {
                    if(a != input_left[l][2 * j] || b != input_right[l][2 * j] || c != input_left[l][2 * j + 1] || d != input_right[l][2 * j + 1]) {
                        cout << "wrong f evaluation" << endl;
                    }
                }
                res = add_modp(res, add_modp(mul_modp(a, b), mul_modp(c, d)));
            }
            if(res != eval_p_poly[l]) {
                cout << "wrong p evaluation, index: " << l << endl;
            }
            // eval_p_poly[l] = res;
        }

        //generate proof
        begin_time = clock();
        vector<uint64_t> ss1(2 * k - 1), ss2(2 * k - 1);
        uint64_t temp;
        for(int i = 0; i < 2 * k - 1; i++) {
            ss1[i] = get_rand();
            if(eval_p_poly[i] > ss1[i]) {
                temp = eval_p_poly[i] - ss1[i];
            }
            else {
                temp = PR - ss1[i] + eval_p_poly[i];
            }
            ss2[i] = temp;
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
        // r = generate_challenge();
        r = rands[cnt];
        cout << "r: " << r << endl; 
        eval_base = evaluate_bases(k, r);
        // for(int j = 0; j < k; j++) {
        //     cout << "evaluate_bases[" << j << "]: " << eval_base[j] << endl;
        // }
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

Proof prove_and_gate(uint64_t _party_id, uint64_t** input_left, uint64_t** input_right, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid, uint64_t* rands) {
    return fliop(input_left, input_right, var, copy, k, sid, rands);
}

struct VerMsg {
    vector<uint64_t> p_eval_ksum_ss;
    vector<uint64_t> p_eval_r_ss;
    vector<uint64_t> input_left_ss;
    vector<uint64_t> input_right_ss;
};

VerMsg gen_vermsg(vector<vector<uint64_t>> p_eval_ss, uint64_t** input_left, uint64_t** input_right, uint64_t** input_mono, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid, uint64_t* rands) {
    uint64_t L = var;
    uint64_t T = copy;
    uint64_t s = (T - 1) / k + 1;
    // uint64_t eta = generate_challenge();
    uint64_t eta = rands[0];
    cout << "in gen_vermsg()" << endl;
    
    // Prepare Input
    begin_time = clock();
    uint64_t eta_power = 1;
    for(int i = 0; i < k; i++) {
        for(int j = 0; j < s; j++) {
            input_left[i][2 * j] = mul_modp(input_left[i][2 * j], eta_power);
            input_left[i][2 * j + 1] = mul_modp(input_left[i][2 * j + 1], eta_power);
            input_mono[i][j] = mul_modp(input_mono[i][j], eta_power);
            eta_power = mul_modp(eta_power, eta);
        }
    }
    finish_time = clock();
    cout<<"Prepare Input Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    uint64_t* eval_base;
    uint64_t r, s0, index, cnt = 0;
    uint128_t temp_result;

    uint64_t len = log(T) / log(k) + 2;
    // uint64_t* p_eval_sum = new uint64_t[len];
    // uint64_t* f_evals_r = new uint64_t[L * len];
    vector<uint64_t> p_eval_ksum_ss(len);
    vector<uint64_t> p_eval_r_ss(len + 1);
    vector<uint64_t> input_left_ss(2 * k);
    vector<uint64_t> input_right_ss(2 * k);

    // Compute share of the monomial's polynomial evaluation at k - 1
    begin_time = clock();
    temp_result = 0;
    for(int i = 0; i < k; i++) {
        temp_result += bacth_sum_modp(input_mono[i], s);
    }
    p_eval_r_ss[0] = modp_128(temp_result);
    finish_time = clock();
    cout<<"Compute Monomial Time = "<<double(finish_time-begin_time)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    uint64_t* eval_base_k = evaluate_bases(k, k - 1);

    s *= 2;
    while(true)
    {
        cout<<"s : "<<s<<endl;
        cout<<"k : "<<k<<endl;

        // Compute share of sum of p's evaluations over [0, k - 1]
        uint64_t res = 0;
        for(int j = 0; j < k; j++) { // Assume k < 64
            res += p_eval_ss[cnt][j];
        }
        p_eval_ksum_ss[cnt] = modp(temp_result);

        r = rands[cnt + 1];
        cout << "r: " << r << endl; 
        eval_base = evaluate_bases(2 * k - 1, r);
        temp_result = 0;
        for(int j = 0; j < 2 * k - 1; j++) {
            temp_result += ((uint128_t) eval_base[j]) * ((uint128_t) p_eval_ss[cnt][j]);
        }
        p_eval_r_ss[cnt + 1] = modp_128(temp_result);

        // Compute share of p's evaluation at r
        begin_time = clock();
        // r = generate_challenge();
        eval_base = evaluate_bases(k, r);
        // for(int j = 0; j < k; j++) {
        //     cout << "evaluate_bases[" << j << "]: " << eval_base[j] << endl;
        // }
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

        if(s == 1) 
        {
            break;
        }

        cnt++;
    }
    cout << "out of loop" << endl;

    for(int i = 0; i < k; i++) {
        input_left_ss[i] = input_left[i][0];
        input_left_ss[2 * i] = input_left[i][1];
        input_right_ss[i] = input_right[i][0];
        input_right_ss[2 * i] = input_right[i][1];
    }

    VerMsg vermsg = {
        p_eval_ksum_ss,
        p_eval_r_ss,
        input_left_ss,
        input_right_ss
    };
    cout << "exiting gen_vermsg()" << endl;
    return vermsg;
}

bool verify_and_gates(vector<vector<uint64_t>> p_eval_ss, uint64_t** input_left, uint64_t** input_right, uint64_t** input_mono, VerMsg other_vermsg, uint64_t var, uint64_t copy, uint64_t k, uint64_t sid, uint64_t* rands) {
    uint64_t L = var;
    uint64_t T = copy;
    uint64_t s = (T - 1) / k + 1;
    uint64_t len = log(T) / log(k) + 1;
    
    VerMsg self_vermsg = gen_vermsg(p_eval_ss, input_left, input_right, input_mono, var, copy, k, sid, rands);
    
    for(int i = 0; i < len; i++) {
        uint64_t p_eval_ksum = add_modp(self_vermsg.p_eval_ksum_ss[i], other_vermsg.p_eval_ksum_ss[i]);
        uint64_t p_eval_r = add_modp(self_vermsg.p_eval_r_ss[i], other_vermsg.p_eval_r_ss[i]);
        if(p_eval_ksum != p_eval_r) {
            cout << i << "-th sum check didn't pass" << endl;
            return false;
        }
    }
    uint64_t* last_input_left = new uint64_t[k];
    uint64_t* last_input_right = new uint64_t[k];
    for(int i = 0; i < k; i++) {
        last_input_left[i] = add_modp(self_vermsg.input_left_ss[i], other_vermsg.input_left_ss[i]);
        last_input_left[2 * i] = add_modp(self_vermsg.input_left_ss[2 * i], other_vermsg.input_left_ss[2 * i]);
        last_input_right[i] = add_modp(self_vermsg.input_right_ss[i], other_vermsg.input_right_ss[i]);
        last_input_right[2 * i] = add_modp(self_vermsg.input_right_ss[2 * i], other_vermsg.input_right_ss[2 * i]);
    }
    uint64_t res = inner_productp(last_input_left, last_input_right, 2 * k);
    uint64_t last_p_eval_r = add_modp(self_vermsg.p_eval_r_ss[len], other_vermsg.p_eval_r_ss[len]);
    if(res != last_p_eval_r) {
        cout << "last check didn't pass" << endl;
        return false;
    }
    
    return true;
}

void shape(uint64_t** input, uint64_t L, uint64_t T, uint64_t k, uint64_t** &input_left, uint64_t** &input_right, uint64_t** &input_mono) {
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

    uint64_t* meta_mono = new uint64_t[s * k];
    input_mono = new uint64_t*[k];
    for(int i = 0; i < k; i++) {
        input_mono[i] = meta_mono + i * s;
        for(int j = 0; j < s; j++) {
            if(i * s + j >= T) {
                input_mono[i][j] = 0;
            }
            else {
                input_mono[i][j] = input[4][i * s + j];
            }
        }
    }
}

int main() {

    // if (!test_evaluate_bases()) {
    //     return 0;
    // }

    uint64_t T = 100;
    uint64_t L = 5;
    uint64_t k = 4;
    uint64_t _party_id = 1;
    srand((unsigned)time(NULL));
    uint64_t sid = get_rand();

    uint64_t cnt = log(T)/log(k) + 2;
    cout<<"cnt : "<< cnt << endl;
    uint64_t* rands = new uint64_t[cnt];
    for(int i = 0; i < cnt; i++) {
        rands[i] = get_rand();
    }

    // Generating satisfying inputs
    uint64_t** input = new uint64_t*[L];
    for(int i = 0; i < L - 1; i++) {
        input[i] = new uint64_t[T];
        for(int j = 0; j < T; j++) {
            input[i][j] = get_rand();
        }
    }
    input[L - 1] = new uint64_t[T];
    for(int j = 0; j < T; j++) {
        uint128_t temp_res = (uint128_t)input[0][j] * (uint128_t)input[1][j] + (uint128_t)input[2][j] * (uint128_t)input[3][j];   
        input[L - 1][j] = modp_128(temp_res);
    }

    // Distributing shares of inputs
    uint64_t** input_ss1 = new uint64_t*[L];
    for(int i = 0; i < L; i++) {
        input_ss1[i] = new uint64_t[T];
        for(int j = 0; j < T; j++) {
            input_ss1[i][j] = get_rand();
        }
    }

    uint64_t** input_ss2 = new uint64_t*[L];
    for(int i = 0; i < L; i++) {
        input_ss2[i] = new uint64_t[T];
        for(int j = 0; j < T; j++) {
            input_ss2[i][j] = sub_modp(input[i][j], input_ss1[i][j]);
        }
    }

    cout << "Checking inputs..." << endl;
    // Check inputs
    for(int j = 0; j < T; j++) {
        uint64_t a, b, c, d, e;
        a = add_modp(input_ss1[0][j], input_ss2[0][j]);
        b = add_modp(input_ss1[1][j], input_ss2[1][j]);
        c = add_modp(input_ss1[2][j], input_ss2[2][j]);
        d = add_modp(input_ss1[3][j], input_ss2[3][j]);
        e = add_modp(input_ss1[4][j], input_ss2[4][j]);
        uint64_t res = add_modp(mul_modp(a, b), mul_modp(c, d));
        if(res != e) {
            cout << "unsatisfied inputs" << endl;
            return 0;
        }
    }
    cout << "inputs satisfied" << endl;


    uint64_t** input_left, **input_right, **input_mono;
    uint64_t** input_left_ss1, **input_right_ss1, **input_mono_ss1;
    uint64_t** input_left_ss2, **input_right_ss2, **input_mono_ss2;

    shape(input, L, T, k, input_left, input_right, input_mono);
    shape(input_ss1, L, T, k, input_left_ss1, input_right_ss1, input_mono_ss1);
    shape(input_ss2, L, T, k, input_left_ss2, input_right_ss2, input_mono_ss2);
    
    cout << "Checking inputs after shape()..." << endl;
    // Check inputs
    for(int i = 0; i < k; i++) {
        uint64_t s = (T - 1) / k + 1;
        uint64_t a, b, c, d, e;
        for(int j = 0; j < s; j++) {
            a = add_modp(input_left_ss1[i][2 * j], input_left_ss2[i][2 * j]);
            b = add_modp(input_right_ss1[i][2 * j], input_right_ss2[i][2 * j]);
            c = add_modp(input_left_ss1[i][2 * j + 1], input_left_ss2[i][2 * j + 1]);
            d = add_modp(input_right_ss1[i][2 * j + 1], input_right_ss2[i][2 * j + 1]);
            e = add_modp(input_mono_ss1[i][j], input_mono_ss2[i][j]);
            uint64_t res = add_modp(mul_modp(a, b), mul_modp(c, d));
            if(res != e) {
                cout << "unsatisfied inputs" << endl;
                return 0;
            }
        }
    }
    cout << "inputs satisfied after shape()" << endl;

    cout<<"T: "<<T<<endl;
    clock_t start, end;

    start = clock();
    Proof proof = prove_and_gate(_party_id, input_left, input_right, L, T, k, sid, rands);
    end = clock();
    cout<<"Total Proving Time = "<<double(end-start)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;

    start = clock();
    VerMsg other_vermsg = gen_vermsg(proof.p_coeffs_ss1, input_left_ss1, input_right_ss1, input_mono_ss1, L, T, k, sid, rands);
    bool res = verify_and_gates(proof.p_coeffs_ss2, input_left_ss2, input_right_ss2, input_mono_ss2, other_vermsg, L, T, k, sid, rands);
    end = clock();
    cout<<"Total Verification Time = "<<double(end-start)/CLOCKS_PER_SEC * 1000<<"ms"<<endl;
    cout<<"Verified = "<<res<<endl;

    return 0;
}

