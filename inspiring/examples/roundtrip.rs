//! End-to-end roundtrip example. The Phase 4 scaffold is intentionally
//! tiny: it imports the public API surface so a `cargo build --example
//! roundtrip` exercises every re-export in `lib.rs`.
//!
//! Phase 9 will replace this file with a worked example that
//!
//! 1. constructs an [`inspiring::RlweParams`] for paper Table 5 set 1,
//! 2. encrypts `d = 1024` random `m_k ∈ Z_p` values as LWE ciphertexts,
//! 3. preprocesses `(A, K_g, K_h)` once,
//! 4. calls [`inspiring::pack`] and decrypts under `s̃`,
//! 5. asserts `|decrypted − Σ m_k X^k|_∞ < Δ/2`.

use inspiring::{GadgetParams, InspiringError};

fn main() -> Result<(), InspiringError> {
    let _gadget = GadgetParams {
        bits_per: 4,
        ell: 8,
    };

    println!("inspiring crate skeleton is wired up. Phase 5+ implementation pending.");
    Ok(())
}
