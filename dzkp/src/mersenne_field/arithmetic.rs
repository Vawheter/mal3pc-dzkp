#![allow(non_snake_case)]

use rand::Rng;

static MERSENNE_PRIME_EXP: usize = 61;
static SHIFT: usize = 3;
static PR: u64 = 2305843009213693951; // 2^61 - 1

pub fn modp(a: u64) -> u64 {
    let i: u64 = (a & PR) + (a >> MERSENNE_PRIME_EXP);
	if i >= PR { i - PR } else { i }
}

// pub fn neg_modp(a: u64) -> u64 {
//     if a > PR {
//         PR - (a - PR)
//     } else {
//         PR - a
//     } 
// }

pub fn neg_modp(a: u64) -> u64 {
    PR - a
}

pub fn add_modp(a: u64, b: u64) -> u64 {
    let res: u64 = a + b; // both less than 2^61 - 1 , thus won't overflow 2^64
    if res >= PR { res - PR } else { res }
}

pub fn sub_modp(a: u64, b: u64) -> u64 {
    if a > b {
        a - b
    }
    else {
        PR - (b - a)
    }
}

// pub fn mul_modp(a: u64, b: u64) -> u64 {
//     let cc: u128 = a as u128 * b as u128;
//     let c = (cc >> 64) as u64;
// 	   let e = cc as u64;
//     let res: u64 = (e & PR) + ((e >> MERSENNE_PRIME_EXP) ^ (c << SHIFT));
// 	   if res >= PR { res - PR } else { res }
// }

pub fn mul_modp(a: u64, b: u64) -> u64 {
    let cc: u128 = a as u128 * b as u128;
    let res = ((cc >> 61) + (cc & PR as u128)) as u64;
    if res >= PR { res - PR } else {res}
}

fn egcd(a: u64, b: u64) -> (u64, u64, u64) {
    if a > b {
        return egcd(b, a);
    }
    if a == 0 {
        (b, 0, 1)
    } else {
        let (g, x, y) = egcd(b % a, a);
        // g = (b % a ) * x + y * a
        // b % a = b - b/a * a
        // let r = (b / a) * x;
        let t: u64 = b / a;
        // assert_eq!(b % a, b - t * a);
        let r = match u64::overflowing_mul(t, x) {
            (v, false) => v,
            (_, true) => mul_modp(t, x),
        };
        if r < y {
            (g, y - r, x)
        } else {
            (g, y + PR - modp(r), x)
        }
    }
}

/// Compute the inverse by means of the Extended Eucleadean Algorithm
pub fn inverse(a: u64) -> u64 {
    let (g, x, _) = egcd(a, PR);
    assert_eq!(g, 1);
    modp(x)
}

pub fn rand_modp<R: Rng>(rng: &mut R) -> u64 {
    rng.gen_range(0, PR)
}

pub fn inner_productp(input_left: &Vec<u64>, input_right: &Vec<u64>) -> u64 {
    let len_left = input_left.len();
    let len_right = input_right.len();
    let bound = (1 << (128 - 2 * MERSENNE_PRIME_EXP)) - 1;
    assert_eq!(bound, 63);
    assert_eq!(len_left, len_right);
    let mut i = 0;
    let mut result : u128 = 0;
    loop{
        let start = i;
        let end = if i + bound < len_left {i + bound} else {len_left};
        for j in start..end {
            result += input_left[j] as u128 * input_right[j] as u128;
        }
        let higher = (result >> (2 * MERSENNE_PRIME_EXP)) as u64;
        let middle = ((result >> MERSENNE_PRIME_EXP) as u64 & PR) as u64;
        let lower = (result as u64 & PR) as u64;
        result = modp(higher + middle + lower) as u128;
        i = end;
        if i == len_left {
            break;
        }
    }
    result as u64
}