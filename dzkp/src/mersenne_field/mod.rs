#![allow(non_snake_case)]

pub mod arithmetic;
pub use arithmetic::*;

#[cfg(test)]
mod tests {
    #[test]
    fn test_modp() {
        use crate::mersenne_field::modp;
        use rand::{thread_rng, Rng};
        
        let PR: u64 = 2305843009213693951; // 2^61 - 1
        let mut rng = thread_rng();
        let n = 100;
        
        for _ in 0..n {
            let a: u64 = rng.gen();
            assert_eq!(a % PR, modp(a));
        }
    }

    #[test]
    fn test_add_sub_modp() {
        use crate::mersenne_field::{rand_modp, add_modp, sub_modp};
        use rand::thread_rng;
        
        let mut rng = thread_rng();
        let n = 100;
        
        for _ in 0..n {
            let a: u64 = rand_modp(&mut rng);
            let b: u64 = rand_modp(&mut rng);
            let c: u64 = add_modp(a, b);
            let d :u64 = sub_modp(c, b);
            let e :u64 = sub_modp(c, a);
            assert_eq!(d, a);
            assert_eq!(e, b);
        }
    }

    #[test]
    fn test_mul_inverse() {
        use crate::mersenne_field::{rand_modp, mul_modp, inverse};
        use rand::thread_rng;
        
        let mut rng = thread_rng();
        let n = 100;
        
        for _ in 0..n {
            let a: u64 = rand_modp(&mut rng);
            let b = inverse(a);
            let c = mul_modp(a, b);
            assert_eq!(c, 1);
        }
    }
}