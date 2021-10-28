import howSlow from "howslow";
import {
  commit,
  proveEval,
  randomFieldElement,
  validateEval,
} from "../src/inner-product-modp.js";

howSlow(
  "inner-product",
  async (start, stop) => {
    let f = Array(1000).fill(0).map(randomFieldElement);

    start("commit");
    let comf = commit(f);
    stop();
    let z = randomFieldElement();

    start("prove");
    let { fz, proof } = await proveEval(f, z);
    stop();

    start("validate");
    let ok = await validateEval(comf, z, fz, proof);
    stop();
  },
  { numberOfRuns: 0, numberOfWarmups: 0 }
);

// BigInt.prototype.toJSON = function () {
//   return Buffer.from(bigIntArrayToBytes([this], basis.byteLength)).toString(
//     "base64"
//   );
// };
// console.log("proof length:", JSON.stringify([comf, proof, z, fz]).length);
// console.log("poly length:", JSON.stringify(f).length);
