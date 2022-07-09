#![allow(non_snake_case)]

pub mod poly;
pub mod interpolation;

pub use poly::*;
pub use interpolation::*;

#[cfg(test)]
mod tests {
    #[test]
    fn test_lag_bases() {
        use crate::polynomial::get_lagrange_bases;

        let n = 100u64;
        let lag_bases = get_lagrange_bases(n);
        for i in 0..n {
            let y = lag_bases[i as usize].evaluate(i);
            assert_eq!(y, 1);
            for j in 0..n {
                if j != i {
                    let y = lag_bases[i as usize].evaluate(j);
                    assert_eq!(y, 0);
                }
            }
        }
    }

    #[test]
    fn test_interpolation() {
        use crate::mersenne_field::rand_modp;
        use crate::polynomial::*;
        use rand::thread_rng;

        let mut rng = thread_rng();
        let n = 100u64;
        
        let evals: Vec<u64> = (0..n).map(|_| rand_modp(&mut rng)).collect();
        let lag_bases = get_lagrange_bases(n);
        let poly = interpolate(&lag_bases, &evals);

        for i in 0..n {
            let y = poly.evaluate(i);
            assert_eq!(evals[i as usize], y);
        }
    }

    #[test]
    fn bench_interpolation() {
        use crate::mersenne_field::*;
        use crate::polynomial::*;
        use rand::thread_rng;
        use std::time::Instant;

        let mut rng = thread_rng();
        let T: usize = 100000;
        let k: usize = 16;
        
        let evals: Vec<Vec<u64>> = (0..T).map(|_| (0..k).map(|_| rand_modp(&mut rng)).collect()).collect();
        let lag_bases = get_lagrange_bases(k as u64);

        // let mut polys: Vec<Poly> = Vec::new();
        let start = Instant::now();
        for i in 0..T {
            let poly = interpolate(&lag_bases, &evals[i]);
            // polys.push(poly);
        }
        let time1 = start.elapsed();
        // for i in 0..T {
        //     println!("polys: {:?}", polys[i].coeffs.len());
        // }

        let evals2: Vec<Vec<u64>> = (0..T).map(|_| (0..k).map(|_| rand_modp(&mut rng)).collect()).collect();
        let mut res = 0u64;
        let start = Instant::now();
        for i in 0..T {
            for j in 0..k {
                res = add_modp(res, mul_modp(evals[i][j], evals2[i][j]));
            }
        }
        let time2 = start.elapsed();
        println!("res: {:?}", res);

        println!("T: {:?}", T);
        println!("k: {:?}", k);
        println!("Interpolate time: {:?}", time1);
        println!("Dot prod time: {:?}", time2);
    }
}