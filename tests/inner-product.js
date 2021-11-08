import test from "ava";
import { P, randomCurvePoint, toBigInts, toPoint } from "../src/bls12-381.js";
import basis from "../src/basis-256-10.js";
import {
  commit,
  proveEval,
  scalarProdGroup,
  validateEval,
} from "../src/inner-product.js";

test("valid", async (t) => {
  let f = [10n, -5n, 213n, 0n, 87691n, 1n];

  let comf = commit(f);
  let z = 11398476n;

  let { fz, proof } = await proveEval(f, z);

  let ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === true, "valid proof is validated");

  ok = await validateEval(comf, z, fz + 1n, proof);
  t.assert(ok === false, "proof for wrong f(z) is not validated");

  let actualValue = proof.transcript[2];
  proof.transcript[2] = actualValue + 1n;
  ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === false, "tampered proof is not validated");

  proof.transcript[2] = actualValue;
  proof.a = proof.a + 123n;
  ok = await validateEval(comf, z, fz, proof);
  t.assert(ok === false, "tampered proof is not validated");
});

test("curve", ({ assert, is }) => {
  // random x coordinate in F_p
  let point = randomCurvePoint();

  assert(point.isOnCurve());
  assert(point.isTorsionFree());

  let { x, y } = toBigInts(point);
  is(x, x % P);
  is(typeof y, "bigint");

  point = toPoint(basis.G[0]);
  assert(point.isOnCurve());
  assert(point.isTorsionFree());

  let a = [1n, 2n, 3n, 4n, 5n, 6n, 7n, 8n, 9n, 10n];
  let G = basis.G.slice(0, 10).map(toPoint);
  let aG = scalarProdGroup(a, G);
  assert(!aG.isZero());
  assert(aG.isOnCurve());
  assert(aG.isTorsionFree());
});
