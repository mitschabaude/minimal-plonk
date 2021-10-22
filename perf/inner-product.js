import howSlow from "howslow";
import { commit, proveEval, validateEval } from "../inner-product-modp.js";
// import { bigIntArrayToBytes } from "../bigint.js";

howSlow(
  "inner-product",
  async (start, stop) => {
    let f = [10n, -5n, 213n, 0n, 87691n, 1n];

    start("commit");
    let comf = commit(f);
    stop();
    let z = 11398476n;

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
