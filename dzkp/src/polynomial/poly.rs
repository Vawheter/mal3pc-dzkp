use crate::mersenne_field::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use rand::Rng;

#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct Poly {
    /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
    pub coeffs: Vec<u64>,
}

impl Poly {
    
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    pub fn one() -> Self {
        Self { coeffs: [1u64].to_vec() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|coeff| *coeff == 0u64)
    }

    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self.coeffs.last().map_or(false, |coeff| *coeff != 0u64));
            self.coeffs.len() - 1
        }
    }

    /// Returns the coefficients of `self`
    pub fn coeffs(&self) -> &[u64] {
        &self.coeffs
    }

    fn truncate_leading_zeros(&mut self) {
        while self.coeffs.last().map_or(false, |c| *c == 0u64) {
            self.coeffs.pop();
        }
    }

    pub fn from_coefficients_vec(coeffs: Vec<u64>) -> Self {
        let mut result = Self { coeffs };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();
        // Check that either the coefficients vec is empty or that the last coeff is
        // non-zero.
        assert!(result.coeffs.last().map_or(true, |coeff| *coeff != 0u64));
        result
    }

    pub fn from_coefficients_slice(coeffs: &[u64]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Outputs a univariate polynomial of degree `d` where
    /// each coefficient is sampled uniformly at random.
    pub fn rand<R: Rng>(d: usize, rng: &mut R) -> Self {
        let mut random_coeffs: Vec<u64> = Vec::new();
        for _ in 0..=d {
            random_coeffs.push(rand_modp(rng));
        }
        Self::from_coefficients_vec(random_coeffs)
    }

    fn horner_evaluate(poly_coeffs: &[u64], point: u64) -> u64 {
        poly_coeffs
            .iter()
            .rfold(0u64, move |result, coeff| add_modp(mul_modp(result, point), *coeff))
    }

    #[cfg(not(feature = "parallel"))]
    fn internal_evaluate(&self, point: u64) -> u64 {
        Self::horner_evaluate(&self.coeffs, point)
    }

    #[cfg(feature = "parallel")]
    fn internal_evaluate(&self, point: u64) -> u64 {
        // Horners method - parallel method
        // compute the number of threads we will be using.
        let num_cpus_available = rayon::current_num_threads();
        let num_coeffs = self.coeffs.len();
        let num_elem_per_thread = max(num_coeffs / num_cpus_available, MIN_ELEMENTS_PER_THREAD);

        // run Horners method on each thread as follows:
        // 1) Split up the coefficients across each thread evenly.
        // 2) Do polynomial evaluation via horner's method for the thread's coefficeints
        // 3) Scale the result point^{thread coefficient start index}
        // Then obtain the final polynomial evaluation by summing each threads result.
        let result = self
            .coeffs
            .par_chunks(num_elem_per_thread)
            .enumerate()
            .map(|(i, chunk)| {
                let mut thread_result = Self::horner_evaluate(&chunk, point);
                thread_result *= modp(point.pow(&[(i * num_elem_per_thread) as u64]));
                thread_result
            })
            .sum();
        result
    }

    pub fn evaluate(&self, point: u64) -> u64 {
        if self.is_zero() {
            return 0u64;
        } else if point == 0u64 {
            return self.coeffs[0];
        }
        self.internal_evaluate(point)
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut result = if self.is_zero() {
            other.clone()
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&other.coeffs)
                .for_each(|(a, b)| {
                    // *a += b;
                    *a = add_modp(*a, *b);
                });
            result
        } else {
            let mut result = other.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&self.coeffs)
                .for_each(|(a, b)| {
                    // *a += b;
                    *a = add_modp(*a, *b);
                });
            result
        };
        result.truncate_leading_zeros();
        result
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut result = if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&other.coeffs)
                .for_each(|(a, b)| {
                    // *a -= b;
                    *a = sub_modp(*a, *b);
                });
            result
        } else {
            let mut result = other.clone();
            result
                .coeffs
                .iter_mut()
                .zip(&self.coeffs)
                .for_each(|(a, b)| {
                    // *a += b;
                    *a = sub_modp(*a, *b);
                    *a = neg_modp(*a);
                });
            result
        };
        result.truncate_leading_zeros();
        result
    }

    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            Self::zero()
        } else {
            let mut result = vec![0u64; self.degree() + other.degree() + 1];
            for (i, self_coeff) in self.coeffs.iter().enumerate() {
                for (j, other_coeff) in other.coeffs.iter().enumerate() {
                    // result[i + j] += &(*self_coeff * other_coeff);
                    result[i + j] = add_modp(result[i + j], mul_modp(*self_coeff, *other_coeff));
                }
            }
            Self { coeffs: result }
        }
    }

    pub fn cmul(&self, c: u64) -> Self {
        if self.is_zero() || c == 0u64 {
            Self::zero()
        } else {
            let mut result = vec![0u64; self.degree() + 1];
            for (i, self_coeff) in self.coeffs.iter().enumerate() {
                result[i] = mul_modp(*self_coeff, c);
            }
            Self { coeffs: result }
        }
    }
}

#[cfg(test)]
mod tests {
    // #[test]
    // fn test_lag_bases {

    // }
}