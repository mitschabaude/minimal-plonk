import howSlow from "howslow";
import { plonkProve, plonkVerify } from "../src/plonk.js";
import { nextPower2 } from "../src/polynomials.js";

howSlow(
  "inner-product",
  async (start, stop) => {
    let ql = [0n, 0n, 0n, 1n];
    let qr = [0n, 0n, 0n, 1n];
    let qo = [-1n, -1n, -1n, -1n];
    let qm = [1n, 1n, 1n, 0n];
    let qc = [0n, 0n, 0n, 0n];
    let selectors = [ql, qr, qo, qm, qc];

    let circuitLength = 40;
    let n = nextPower2(circuitLength);
    let permutation = [
      [n, n + 1, n + 2, 2 * n],
      [0, 1, 2, 2 * n + 1],
      [3, n + 3, 2 * n + 3, 2 * n + 2],
    ];
    let circuit = { selectors, permutation, circuitLength };

    let a = [3n, 4n, 5n, 9n];
    let b = [3n, 4n, 5n, 16n];
    let c = [9n, 16n, 25n, 25n];
    let witness = [a, b, c];

    start("prove");
    let snark = await plonkProve(circuit, witness);
    stop();

    start("verify");
    let ok = await plonkVerify(circuit, snark);
    stop();

    console.assert(ok);
  },
  { numberOfRuns: 0, numberOfWarmups: 0 }
);
