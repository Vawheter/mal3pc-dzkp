#![allow(non_snake_case)]
#![allow(unused_assignments)]

use crate::mersenne_field::*;
use crate::polynomial::*;
use merlin::Transcript;
use rand::Rng;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write};

use std::time::Instant;

static PR: u64 = 2305843009213693951; // 2^61 - 1

#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof {
    p_coeffs_ss1: Vec<Vec<u64>>,
    p_coeffs_ss2: Vec<Vec<u64>>,
}

fn get_bases(n: u64) -> Vec<Vec<u64>> {
    assert!(n < PR);
    let mut result : Vec<Vec<u64>> = vec![];
    let end = 2 * n - 1;
    for l in n..end {
        let mut row : Vec<u64> = vec![];
        for i in 0..n {
            let mut entry : u64 = 1;
            for j in 0..n {
                if j != i {
                    let denominator = if i > j { i - j } else { neg_modp(j - i) };
                    let numerator = l - j;
                    entry = mul_modp(entry, mul_modp(inverse(denominator), numerator));
                }
            }
            row.push(entry);
        }
        result.push(row);
    }
    result
}

fn evaluate_base(n: u64, r: u64) -> Vec<u64> {
    let mut result : Vec<u64> = vec![];
    for i in 0..n{
        let mut entry : u64 = 1;
        for j in 0..n {
            if j != i {
                let denominator = if i > j { i - j } else { neg_modp(j - i) };
                let numerator = if r > j { r - j } else { neg_modp(j - r)};
                entry = mul_modp(entry, mul_modp(inverse(denominator), numerator));
            }
        }
        result.push(entry);
    }
    result
}

fn matrix_mult(matrix: &Vec<Vec<u64>>, input: &Vec<u64>) -> Vec<u64> {
    let mut result : Vec<u64> = vec![];
    let len = matrix.len();
    for i in 0..len {
        result.push(inner_productp(&matrix[i], &input));
    }
    result
}

fn fliop<R: Rng>(
    inputs: &Vec<Vec<u64>>,
    ks: &[usize], // compression parameter
    sid: usize, // session id
    rng: &mut R,
) -> Proof {
    let L: usize = inputs.len(); // number of variables
    let T: usize = inputs[0].len(); // number of copies

    let mut k: usize = ks[0];
    let mut s0 = T; // s in the last round
    let mut s = (T-1) / k + 1; // s in the current round

    // Add public information --- May have Security Issue
    let mut buf1 = [0u8; 8];
    let mut buf2 = [0u8; 8];

    let mut transcript1 = Transcript::new(b"DZKP_FLIOP_1");
    let mut transcript2 = Transcript::new(b"DZKP_FLIOP_2");
    transcript1.append_message(b"k0", &(ks[0] as u64).to_be_bytes());
    transcript1.append_message(b"k1", &(ks[1] as u64).to_be_bytes());
    transcript1.append_message(b"sid", &(sid as u64).to_be_bytes());
    transcript2.append_message(b"k0", &(ks[0] as u64).to_be_bytes());
    transcript2.append_message(b"k1", &(ks[1] as u64).to_be_bytes());
    transcript2.append_message(b"sid", &(sid as u64).to_be_bytes());

    transcript1.challenge_bytes(b"eta", &mut buf1);
    let eta = u64::from_be_bytes(buf1) & PR;

    // Compute Input
    let input_start = Instant::now();
    let mut eta_power = 1u64;
    let mut input_left = Vec::new();
    for i in 0..s0 {
        input_left.push(mul_modp(eta_power, inputs[0][i]));
        input_left.push(mul_modp(eta_power, inputs[2][i]));
        eta_power = mul_modp(eta_power, eta);
    }

    let mut input_right = Vec::new();
    for i in 0..s0 {
        input_right.push(inputs[1][i]);
        input_right.push(inputs[3][i]);
    }
    let input_time = input_start.elapsed();
    println!("Input Time: {:?}", input_time);

    let mut input_right = Vec::new();
    let mut input_left= Vec::new();
    // Proof
    let mut p_coeffs_ss1: Vec<Vec<u64>> = Vec::new();
    let mut p_coeffs_ss2: Vec<Vec<u64>> = Vec::new();

    // Compute Bases
    let lag_start = Instant::now();
    let base_0 = get_bases(ks[0] as u64);
    let base_1 = get_bases(ks[1] as u64);
    let lag_time = lag_start.elapsed();
    println!("Lagrange time: {:?}", lag_time);
    
    let mut base_ptr = &base_0;

    loop {
        s0 = input_left.len();
        s = (s0 - 1)/k + 1;
        println!("s0: {:?}", s0);
        println!("s: {:?}", s);
        println!("k: {:?}", k);

        // Interpolation
        let interpolation_start = Instant::now();
        let mut eval_left_polys: Vec<Vec<u64>> = Vec::new();
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let mut result = input_left[cur..end].to_vec();
            let num_of_0 = k - (end - cur);
            for j in 0..num_of_0 {
                result.push(0);
            }
            let mut eval = matrix_mult(base_ptr, &result);
            result.append(&mut eval);
            eval_left_polys.push(result);
        }

        let mut eval_right_polys: Vec<Vec<u64>> = Vec::new();
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let mut result = input_right[cur..end].to_vec();
            let num_of_0 = k - (end - cur);
            for j in 0..num_of_0 {
                result.push(0);
            }
            let mut eval = matrix_mult(base_ptr, &result);
            result.append(&mut eval);
            eval_right_polys.push(result);
        }
        let interpolation_time = interpolation_start.elapsed();
        println!("Interpolation time: {:?}", interpolation_time);

        // Compute p(X)
        let poly_mul_start = Instant::now();
        let mut eval_p_poly : Vec<u64> = Vec::new();
        let end = 2 * k - 1;
        for i in 0..end {
            let mut result : u128 = 0;
            let mut j = 0;
            let bound = 63;
            loop{
                let start = j;
                let end = if j + bound < s {j + bound} else {s};
                for j in start..end {
                    result += eval_left_polys[j][i] as u128 * eval_left_polys[j][i] as u128;
                }
                let higher = (result >> 122) as u64;
                let middle = ((result >> 61) as u64 & PR) as u64;
                let lower = (result as u64 & PR) as u64;
                result = modp(higher + middle + lower) as u128;
                j = end;
                if j == s {
                    break;
                }
            }
            eval_p_poly.push(result as u64)
            //let left_vec = (0..s).map(|j| eval_left_polys[j][i]).collect();
            //let right_vec = (0..s).map(|j| eval_right_polys[j][i]).collect();
            //eval_p_poly.push(inner_productp(&left_vec, &right_vec));
        }
        let poly_mul_time = poly_mul_start.elapsed();
        println!("Poly mul time: {:?}", poly_mul_time);

        // Generate proof
        let ss1: Vec<u64> = (0..end).map(|_| rand_modp(rng)).collect();
        let ss2: Vec<u64> = (0..end).map(|i| sub_modp(eval_p_poly[i], ss1[i])).collect();

        transcript1.append_message(b"ss1", unsafe { ss1.align_to::<u8>().1 });
        transcript2.append_message(b"ss2", unsafe { ss2.align_to::<u8>().1 });

        p_coeffs_ss1.push(ss1);
        p_coeffs_ss2.push(ss2);

        if s == 1 {
            break;
        }

        // Verifiers need to exchange their challenges once
        transcript1.challenge_bytes(b"r1", &mut buf1);
        transcript2.challenge_bytes(b"r2", &mut buf2);
        let r1 = u64::from_be_bytes(buf1) & PR;
        let r2 = u64::from_be_bytes(buf2) & PR;

        let r = add_modp(r1, r2);

        let evaluation_start = Instant::now();
        
        // Compute new inputs
        let eval_r_base = evaluate_base(k as u64, r);
        println!("Len of Base: {:?}", eval_r_base.len());
        let mut new_input_left = vec![];
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let len = end - cur;
            let mut right_vec = input_left[cur..end].to_vec();
            let num_of_0 = k - (end - cur);
            for j in 0..num_of_0 {
                right_vec.push(0);
            }
            new_input_left.push(inner_productp(&eval_r_base, &right_vec));
        }
        input_left = new_input_left;

        let mut new_input_right = vec![];
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let len = end - cur;
            let mut right_vec = input_right[cur..end].to_vec();
            let num_of_0 = k - (end - cur);
            for j in 0..num_of_0 {
                right_vec.push(0);
            }
            new_input_right.push(inner_productp(&eval_r_base, &right_vec));
        }
        input_right = new_input_right;
        let evaluation_time = evaluation_start.elapsed();
        println!("Evaluation time: {:?}", evaluation_time);

        k = ks[1];
        base_ptr = &base_1;
    }

    Proof {
        p_coeffs_ss1,
        p_coeffs_ss2,
    }
}

pub fn prove_and_gates<R: Rng>(
    _party_id: usize,
    inputs: &Vec<Vec<u64>>,
    ks: &[usize], // compression parameter
    sid: usize, // session id
    rng: &mut R,
) -> Proof {
    fliop(inputs, ks, sid, rng) //, &circuit1)
    // match party_id {
    //     0 => fliop(inputs, k, sid, rng, &circuit0),
    //     1 => fliop(inputs, k, sid, rng, &circuit1),
    //     2 => fliop(inputs, k, sid, rng, &circuit2),
    // }
}