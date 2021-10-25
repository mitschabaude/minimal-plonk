import howSlow from "howslow";
import {
  commit,
  proveEval,
  randomFieldElement,
  validateEval,
} from "../inner-product-modp.js";

howSlow(
  "inner-product",
  async (start, stop) => {
    let f = Array(10).fill(0).map(randomFieldElement);

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
  { numberOfRuns: 10 }
);

// BigInt.prototype.toJSON = function () {
//   return Buffer.from(bigIntArrayToBytes([this], basis.byteLength)).toString(
//     "base64"
//   );
// };
// console.log("proof length:", JSON.stringify([comf, proof, z, fz]).length);
// console.log("poly length:", JSON.stringify(f).length);
