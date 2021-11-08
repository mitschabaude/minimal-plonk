# minimal-plonk

This is a minimal implementation of [PLONK](https://eprint.iacr.org/2019/953), instantiated with the inner product argument PCS on BLS12-381. It was created for myself as a learning experience, leaves out many optimizations and protocol features like zero knowledge, and is neither efficient nor secure.

- **Inner product argument:** [./inner-product.js](https://github.com/mitschabaude/zktoys/blob/main/inner-product.js) implements the _polynomial commitment scheme (PCS)_ described [here](https://www.cryptologie.net/article/528/what-is-an-inner-product-argument-part-1/), [here](https://dankradfeist.de/ethereum/2021/07/27/inner-product-arguments.html) and [here](https://doc-internal.dalek.rs/bulletproofs/notes/inner_product_proof/index.html).
