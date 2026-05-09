

# InsPIRing

## Key Insights

### LWE-to-RLWE Zeroeth Coefficient

When packing, represent LWE in RLWE form.

This transforms:
```
(a, b = -<a, s> + m)
```

where m is a scalar

to 

```
(ã, b̃ = -ã · s̃ + m̃)
```

where m̃ is a polynomial of a degree d. The message is in the zeroeth coefficient. The rest of the coefficients hold garbage:

```
m + 0X + 0X^2 + ... + 0X^(d-1)
```

In CDKS, we recursively FFT-merge ciphertexts, applying trace functions at every recursive level to
shift-and-zero garbage coefficients, placing messages in their respective terms.

**In InsPIRing, the goal is to transform the ciphertext so that the message becomes only the constant coefficients**
```
m̂(X) = m
```

### The Trace Trick

For a polynomial,
```
p(X) = c0 + c1X + c2X^2 + ... + c(d-1)X^(d-1)
```

Lemma 1 says that if you sum certain Galois automorphisms of p, all non-constant coefficients cancel out, leaving:
```
Tr(p) = d · c0
```

So the trace acts like a “constant-term extractor,” up to a factor of d.

That is why the construction uses d^-1: it scales things so that after applying the trace idea, you get c0, not d · c0.

### Building the Intermediate Ciphertext

The new cuphertext has the form:
```
(â, b̃)
```

where `â` is not a singly polynomial. It is a vector of d polynomials:
```
â ∈ Rq^d
```

The construction fills `â` using authomorphisms of the original RLWE random component `ã`.

For the first half:
```
â[j] = d^-1 · τg^j(ã)
```

For the second half:
```
â[d/2 + j] = τh(â[j])
```

So â is a structured vector of Galois-transformed copies of ã.

The secret key is constructed in the matching way:
```
ŝ[j] = τg^j(s̃)
ŝ[d/2 + j] = τh(ŝ[j])
```

This gives a secret key vector whose components are highly correlated. They are not independent random secrets. They are all automorphic images of the same original secret.

That correlation is intentional. This structure is needed later, where two key-switching matrices are used.



## Open Questions

- What are the security implications on this new structure? We clearly make the encryption more structured than before.


## Implications

- InsPIRing is strictly worse than CDKS in-terms of pre-processing. The main idea of InsPIRing is to offload more online computation to offline. That may be constraining for databases needing a fast rebuild.
