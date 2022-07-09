#![allow(non_snake_case)]

pub mod prover;
pub mod prover2;

pub use prover::*;
pub use prover2::*;

#[cfg(test)]
mod tests {
    #[test]
    fn bench_prove_and_gates() {
        use crate::mersenne_field::rand_modp;
        use crate::binary_dzkp::prove_and_gates;
        use std::time::Instant;
        use rand::{thread_rng, Rng};
        use ark_serialize::CanonicalSerialize;

        let T: usize = 10000;
        let L: usize = 5;
        let ks = [8, 8]; 
        let party_id: usize = 1;

        println!("T: {}", T);
    
        let mut rng = &mut thread_rng();
        let sid: usize = rng.gen();

        let inputs = (0..L).map(|_| (0..T).map(|_| rand_modp(&mut rng)).collect()).collect();
        // let inputs = (0..L).map(|_| (0..T).map(|_| rng.gen_range(0,2)).collect()).collect();

        let prove_start = Instant::now();
        let proof = prove_and_gates(party_id, &inputs, &ks, sid, &mut rng);
        let prove_time = prove_start.elapsed();
        println!("Proving time: {:?}", prove_time);

        let mut proof_bytes = vec![];
        proof.serialize(&mut proof_bytes).unwrap();
        println!("Proof length: {} bytes", proof_bytes.len());
    }
}