# minimal-plonk

This is a minimal implementation of [PLONK](https://eprint.iacr.org/2019/953) in JavaScript, instantiated with the inner product argument PCS on BLS12-381. It was created for myself as a learning experience, leaves out protocol features like zero knowledge, doesn't provide an actual way to create circuits and witnesses, and is neither efficient nor secure.

## TODOs

The first thing to address would be the awful speed. Currently the runtime of both prover & verifier is completely dominated (>98%) by scalar multiplications. These should be made much more efficient by

- using linearization/batching techniques
- implementing a smaller curve, e.g. [Pasta](https://electriccoin.co/blog/the-pasta-curves-for-halo-2-and-beyond/)
- using multi scalar multiplication (MSM) techniques, e.g. bucket method / Pippenger
- implementing core routines in Wasm (at least double / add operations, better the entire MSM)

After that, it would be interesting to start testing some non-trivial circuits and think about how to easily create them.

## Details

- **Inner product argument:** [./src/inner-product.js](https://github.com/mitschabaude/minimal-plonk/blob/main/src/inner-product.js) implements the _polynomial commitment scheme (PCS)_ described [here](https://www.cryptologie.net/article/528/what-is-an-inner-product-argument-part-1/), [here](https://dankradfeist.de/ethereum/2021/07/27/inner-product-arguments.html) and [here](https://doc-internal.dalek.rs/bulletproofs/notes/inner_product_proof/index.html).

- **PLONK:** [./src/plonk.js](https://github.com/mitschabaude/minimal-plonk/blob/main/src/plonk.js) implements the Prover / Verifier algorithms of PLONK over IPA, with no special gates, no plookup, no zero knowledge, very inefficient use of the PCS.

## Test

To run tests in watch mode:

```sh
npx ava --watch
```

To inspect performance in Chrome:

```
npx chrode perf/plonk.js --no-headless
```
