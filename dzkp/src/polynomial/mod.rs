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
        use crate::mersenne_field::modp;
        use crate::polynomial::*;
        use rand::{thread_rng, Rng};

        let mut rng = thread_rng();
        let n = 100u64;
        
        let evals: Vec<u64> = (0..n).map(|_| modp(rng.gen())).collect();
        let lag_bases = get_lagrange_bases(n);
        let poly = interpolate(&lag_bases, &evals);

        for i in 0..n {
            println!("i: {}", i);
            let y = poly.evaluate(i);
            assert_eq!(evals[i as usize], y);
        }
    }
}