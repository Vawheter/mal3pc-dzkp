use crate::mersenne_field::*;
use crate::polynomial::Poly;

static PR: u64 = 2305843009213693951; // 2^61 - 1

pub fn get_lagrange_bases(n: u64) -> Vec<Poly> {
    assert!(n < PR);
    let mut lang_polys = Vec::new();
    for i in 0..n {
        let mut li_poly = Poly::from_coefficients_vec([1u64].to_vec());
        // let range = if i == 0 { (1..n) } else { (0..i - 1).chain(i + 1..n) };
        for j in 0..n {
            if j != i {
                let delta = if i > j { i - j } else { PR - (j - i) };
                let c = inverse(delta);
                let c0 = mul_modp(c, PR - j);
                let pij = Poly::from_coefficients_vec([c0, c].to_vec());
                li_poly = li_poly.mul(&pij);
            }
        }
        lang_polys.push(li_poly);
    }
    lang_polys
}

pub fn interpolate(lag_bases: &Vec<Poly>, evals: &Vec<u64>) -> Poly {
    let n = evals.len();
    assert_eq!(lag_bases.len(), n);
    let mut coeffs = vec![0u64; n];
    for i in 0..n {
        for j in 0..n {
            coeffs[j] = add_modp(coeffs[j], mul_modp(lag_bases[i].coeffs[j], evals[i]));
        }
    }
    Poly::from_coefficients_vec(coeffs)
}