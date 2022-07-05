pub mod poly;
pub mod lagrange;

pub use poly::*;
pub use lagrange::*;

#[cfg(test)]
mod tests {
    #[test]
    fn test_lag_bases() {
        use crate::polynomial::{Poly, get_lagrange_bases};

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
}