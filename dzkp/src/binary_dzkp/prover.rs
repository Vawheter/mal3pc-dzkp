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

fn fliop<R: Rng>(
    inputs: &Vec<Vec<u64>>,
    ks: &[usize], // compression parameter
    sid: usize, // session id
    rng: &mut R,
    // circuit: &dyn Fn(&Vec<Vec<Poly>>, &Vec<u64>) -> Poly,
) -> Proof {
    let L: usize = inputs.len(); // number of variables
    let T: usize = inputs[0].len(); // number of copies

    let mut k: usize = ks[0];

    let S: usize = (T as f64 / k as f64).ceil() as usize; // rounding down
    let mut s0 = T; // s in the last round
    let mut s = S; // s in the current round

    // Inputs to circuits in iterations
    let mut f_evals_r = inputs.clone();

    let mut buf1 = [0u8; 8];
    let mut buf2 = [0u8; 8];
    let map = u64::MAX / PR;

    let mut transcript1 = Transcript::new(b"DZKP_FLIOP_1");
    let mut transcript2 = Transcript::new(b"DZKP_FLIOP_2");
    // Add public information
    transcript1.append_message(b"k0", &(ks[0] as u64).to_be_bytes());
    transcript1.append_message(b"k1", &(ks[1] as u64).to_be_bytes());
    transcript1.append_message(b"sid", &(sid as u64).to_be_bytes());
    transcript2.append_message(b"k0", &(ks[0] as u64).to_be_bytes());
    transcript2.append_message(b"k1", &(ks[1] as u64).to_be_bytes());
    transcript2.append_message(b"sid", &(sid as u64).to_be_bytes());

    transcript1.challenge_bytes(b"eta", &mut buf1);
    let eta = u64::from_be_bytes(buf1) / map;

    let mut eta_power = 1u64;
    let mut eta_powers: Vec<u64> = Vec::new();
    for _ in 0..s {
        eta_powers.push(eta_power);
        eta_power = mul_modp(eta_power, eta);
    }

    let mut f_eta_evals_r = eta_powers.clone();
    let mut f_eta_evals_r0 = Vec::new();

    let mut p_coeffs_ss1: Vec<Vec<u64>> = Vec::new();
    let mut p_coeffs_ss2: Vec<Vec<u64>> = Vec::new();

    let lag_start = Instant::now();
    let lag_bases_0 = get_lagrange_bases(ks[0] as u64);
    let lag_bases_1 = get_lagrange_bases(ks[1] as u64);
    let lag_time = lag_start.elapsed();
    println!("Lagrange time: {:?}", lag_time);
    
    let mut lag_ptr = &lag_bases_0;

    loop {
        println!("s: {:?}", s);
        println!("k: {:?}", k);

        // Interpolation
        let interpolation_start = Instant::now();
        let mut f_polys: Vec<Vec<Poly>> = Vec::new();
        for l in 0..L - 1 {
            let mut f_l_polys: Vec<Poly> = Vec::new();
            for i in 0..s {
                let cur = i * k;
                let end = if i == s - 1 { s0 } else { cur + k }; 
                let poly = interpolate(lag_ptr, &f_evals_r[l][cur..end].to_vec());
                f_l_polys.push(poly);
            }
            f_polys.push(f_l_polys);
        }
        let mut agg_input = vec![0u64; k];
        let t = f_evals_r[L - 1].len();
        for i in 0..s {
            let len = if i == s - 1 { t - (s - 1) * k } else { k }; 
            for j in 0..len {
                agg_input[j] = add_modp(agg_input[j], mul_modp(f_evals_r[L - 1][i * k + j], f_eta_evals_r[i]));
            }
        }
        let poly = interpolate(lag_ptr, &agg_input);
        f_polys.push([poly].to_vec());
        let interpolation_time = interpolation_start.elapsed();
        println!("Interpolation time: {:?}", interpolation_time);

        // Compute p(X)
        let poly_mul_start = Instant::now();
        let mut p_poly = Poly::zero();
        for i in 0..s {
            let tmp1 = f_polys[0][i].mul(&f_polys[1][i]);
            let tmp2 = f_polys[2][i].mul(&f_polys[3][i]);
            let res = tmp1.add(&tmp2).cmul(f_eta_evals_r[i]);
            p_poly = p_poly.add(&res);
        }
        p_poly = p_poly.add(&f_polys[L - 1][0]);
        let poly_mul_time = poly_mul_start.elapsed();
        println!("Poly mul time: {:?}", poly_mul_time);

        // Generate proof
        let ss1: Vec<u64> = (0..k).map(|_| rand_modp(rng)).collect();
        let ss2: Vec<u64> = (0..k).map(|i| sub_modp(p_poly.coeffs[i], ss1[i])).collect();

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
        let r1 = u64::from_be_bytes(buf1) / map;
        let r2 = u64::from_be_bytes(buf2) / map;

        let r = add_modp(r1, r2);

        let evaluation_start = Instant::now();
        
        // Compute new inputs
        let mut var_evals_r = Vec::new();
        let lag_evals_r: Vec<u64> = (0..k).map(|i| lag_ptr[i].evaluate(r)).collect();
        // println!("f_evals_r[L - 1].len(): {}", f_evals_r[L - 1].len()); 
        for i in 0..s {
            let len = if i == s - 1 { t - (s - 1) * k } else { k };
            // println!("len: {}", len); 
            // println!("i * k: {}", i * k); 
            var_evals_r.push((0..len).map(|j|
                mul_modp(f_evals_r[L - 1][i * k + j], lag_evals_r[j])
            ).sum());
        }

        f_evals_r = Vec::new();
        for l in 0..L - 1 {
            f_evals_r.push((0..s).map(|i| 
                f_polys[l][i].evaluate(r)
            ).collect());
        }
        f_evals_r.push(var_evals_r);

        // the second round till the last round
        k = ks[1];
        lag_ptr = &lag_bases_1;
        s0 = s;
        s = (s0 as f64 / k as f64).ceil() as usize;

        f_eta_evals_r0 = f_eta_evals_r.clone();
        f_eta_evals_r = vec![];
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let poly = interpolate(&lag_ptr, &f_eta_evals_r0[cur..end].to_vec());
            f_eta_evals_r.push(poly.evaluate(r));
        }
        let evaluation_time = evaluation_start.elapsed();
        println!("Evaluation time: {:?}", evaluation_time);

    }

    Proof {
        p_coeffs_ss1,
        p_coeffs_ss2,
    }
}

/// P0's local computaion
// fn circuit0(input_polys: &Vec<Vec<Poly>>, etas: &Vec<u64>) -> Poly {

// }

/// P1's local computaion
// v0 = \delta_x\delta_y\delta_{z_1}, v1 = r_{x_1}
// v2 = \delta_x, v3 = r_{y_1}r_{xy,1}r_{z_1}
// v4 = \delta_x\delta_y, v5 = v1 = r_{x_1}
// v5 = \delta_y\delta_{z_1}, v6 = v1 = r_{x_1}
// v6 = v2 = \delta_x, v6 = r_{y_1}r_{xy,1}
// v7 = v2 = \delta_x, v7 = r_{y_1}r_{z_1}
// v8 = \delta_y, v9 = v1 = r_{x_1}
// v9 = v2 = \delta_x, v9 = r_{y_1}
// v10 = \delta_x\delta_y\delta_{z_1} + \delta_x\delta_y + \delta_{z_1}
// v11 = r_{xy,1}r_{z_1} - r_{xy,1} - r_{z_1}
// 4(v0v1 - v2v3) - 2(v4v1 + v5v1 + v2v6 + v2v7) + v8v1 + v2v9 + v10 + v11
// (4v0 - 2(v4 + v5) + v8)v1 - (4v3 + 2v6 + 2v7 - v9)v2 + v10 + v11
// v0 = 4v0 - 2(v4 + v5) + v8, v1 = v1
// v2 = 4v3 + 2v6 + 2v7 - v9, v3 = v2
// v4 = v10 + v11
// v0 * v1 + v2 * v3 + v4 = 0


/// P2's local computaion
// fn circuit2(input_polys: &Vec<Vec<Poly>>, etas: &Vec<u64>) -> Poly {

// }

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