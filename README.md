# ZK toys

This is just me trying to implement SNARK algorithms for my own understanding. The code is simple & educational rather than efficient | secure.

What's here so far:

- **Inner product argument PCS:** [./inner-product-modp.js](https://github.com/mitschabaude/zktoys/blob/main/inner-product-modp.js) implements the Bulletproof-style _polynomial commitment scheme (PCS)_ described [here](https://www.cryptologie.net/article/528/what-is-an-inner-product-argument-part-1/), [here](https://dankradfeist.de/ethereum/2021/07/27/inner-product-arguments.html) and [here](https://doc-internal.dalek.rs/bulletproofs/notes/inner_product_proof/index.html).
  - A difference to those descriptions is that I don't use an elliptic curve as the group but the multiplicative group of integers mod p, which is fairly simple to implement from scratch thanks to built-in JS `BigInt`.
  - There's also code to generate a random large prime and random "basis" elements in Z_p.
