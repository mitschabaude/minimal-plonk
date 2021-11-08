import test from "ava";
import { commit, proveEval, validateEval } from "../src/inner-product-modp.js";
// import { bigIntArrayToBytes } from "../bigint.js";

test("valid", async (t) => {
  let f = [10n, -5n, 213n, 0n, 87691n, 1n];

  let comf = commit(f);
  let z = 11398476n;

  let { fz, proof } = await proveEval(f, z);

  let ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === true, "valid proof is validated");

  ok = await validateEval(comf, z, fz + 1n, proof);
  t.assert(ok === false, "proof for wrong f(z) is not validated");

  let actualValue = proof.transcript[0];
  proof.transcript[0] = actualValue + 1n;
  ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === false, "tampered proof is not validated");

  proof.transcript[0] = actualValue;
  proof.a = proof.a + 123n;
  ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === false, "tampered proof is not validated");
});
