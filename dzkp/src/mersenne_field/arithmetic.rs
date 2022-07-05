static MERSENNE_PRIME_EXP: usize = 61;
static SHIFT: usize = 3;
pub static PR: u64 = 2305843009213693951; // 2^61 - 1
static PR_128: u128 = 2305843009213693951; // 2^61 - 1

pub fn modp(a: u64) -> u64 {
    let i: u64 = (a & PR) + (a >> MERSENNE_PRIME_EXP);
	if i >= PR {i - PR} else {i}
}

pub fn add_modp(a: u64, b: u64) -> u64 {
    let res: u64 = a + b; // all less than 2^61 - 1 , thus won't overflow 2^64
    if res >= PR {res - PR} else {res}
}

pub fn sub_modp(a: u64, b: u64) -> u64 {
    if a > b {
        a - b
    }
    else {
        PR - (b - a)
    }
}

pub fn mul_modp(a: u64, b: u64) -> u64 {
    let cc: u128 = a as u128 * b as u128;
    let c = (cc >> 64) as u64;
	let e = cc as u64;
    let res: u64 = (e & PR) + ((e >> MERSENNE_PRIME_EXP) ^ (c << SHIFT));
	if res >= PR {res - PR} else {res}
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