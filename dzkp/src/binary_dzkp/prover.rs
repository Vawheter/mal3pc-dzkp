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
    k: usize, // compression parameter
    sid: usize, // session id
    rng: &mut R,
    circuit: &dyn Fn(&Vec<Vec<Poly>>, &Vec<u64>) -> Poly,
) -> Proof {
    let L: usize = inputs.len(); // number of variables
    let T: usize = inputs[0].len(); // number of copies

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
    transcript1.append_message(b"k", &(k as u64).to_be_bytes());
    transcript1.append_message(b"sid", &(sid as u64).to_be_bytes());
    transcript2.append_message(b"k", &(k as u64).to_be_bytes());
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
    let mut f_eta_polys: Vec<Poly> = Vec::new();

    let mut p_coeffs_ss1: Vec<Vec<u64>> = Vec::new();
    let mut p_coeffs_ss2: Vec<Vec<u64>> = Vec::new();

    let lag_start = Instant::now();
    let lag_bases = get_lagrange_bases(k as u64);
    let lag_time = lag_start.elapsed();
    println!("Lagrange time: {:?}", lag_time);

    loop {
        println!("s: {:?}", s);

        // Interpolation
        let interpolation_start = Instant::now();
        let mut f_polys: Vec<Vec<Poly>> = Vec::new();
        for l in 0..L {
            let mut f_l_polys: Vec<Poly> = Vec::new();
            for i in 0..s {
                let cur = i * k;
                let end = if i == s - 1 { s0 } else { cur + k }; 
                let poly = interpolate(&lag_bases, &f_evals_r[l][cur..end].to_vec());
                f_l_polys.push(poly);
            }
            f_polys.push(f_l_polys);
        }
        let interpolation_time = interpolation_start.elapsed();
        println!("Interpolation time: {:?}", interpolation_time);

        // Compute p(X)
        let poly_mul_start = Instant::now();
        let p_poly = circuit(&f_polys, &f_eta_evals_r);
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
        f_evals_r = Vec::new();
        for l in 0..L {
            f_evals_r.push((0..s).map(|i| 
                f_polys[l][i].evaluate(r)
            ).collect());
        }

        s0 = s;
        s = (s0 as f64 / k as f64).ceil() as usize;

        f_eta_evals_r0 = f_eta_evals_r.clone();
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 { s0 } else { cur + k };
            let poly = interpolate(&lag_bases, &f_eta_evals_r0[cur..end].to_vec());
            f_eta_polys.push(poly);
        }

        f_eta_evals_r = vec![];
        for i in 0..s {
            f_eta_evals_r.push(f_eta_polys[i].evaluate(r));
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
fn circuit1(input_polys: &Vec<Vec<Poly>>, etas: &Vec<u64>) -> Poly
{
    assert_eq!(input_polys.len(), 5); // 5 variables in total
    let s = input_polys[0].len();
    let mut poly = Poly::zero();
    for i in 0..s {
        let tmp1 = input_polys[0][i].mul(&input_polys[1][i]);
        let tmp2 = input_polys[2][i].mul(&input_polys[3][i]);
        let res = tmp1.add(&tmp2).add(&input_polys[4][i]);
        poly = poly.add(&res.cmul(etas[i]));
    }
    poly
}

/// P2's local computaion
// fn circuit2(input_polys: &Vec<Vec<Poly>>, etas: &Vec<u64>) -> Poly {

// }

pub fn prove_and_gates<R: Rng>(
    _party_id: usize,
    inputs: &Vec<Vec<u64>>,
    k: usize, // compression parameter
    sid: usize, // session id
    rng: &mut R,
) -> Proof {
    fliop(inputs, k, sid, rng, &circuit1)
    // match party_id {
    //     0 => fliop(inputs, k, sid, rng, &circuit0),
    //     1 => fliop(inputs, k, sid, rng, &circuit1),
    //     2 => fliop(inputs, k, sid, rng, &circuit2),
    // }
}

