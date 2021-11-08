import bls, { PointG1, Fp } from "@noble/bls12-381";
import { mod } from "./modular-arithmetic.js";
import { randomFieldElement } from "./random-primes.js";
const P = bls.CURVE.P;

export { P, randomCurvePoint, toPoint, toBigInts };

function toPoint({ x, y }) {
  return new PointG1(new Fp(x), new Fp(y));
}
function toBigInts(point) {
  let [{ value: x }, { value: y }] = point.toAffine();
  return { x, y };
}

function randomCurvePoint() {
  let x = randomFieldElement(P, 48);
  let y;
  while (true) {
    // compute y^2 = x^3 + 4
    let ysquare = mod(x ** 3n + 4n, P);

    // try computing square root to get y (works half the time, because half the field elements are squares)
    y = new Fp(ysquare).sqrt();
    if (y !== undefined) {
      break;
    } else {
      // if it didn't work, increase x by 1 and try again
      x = mod(x + 1n, P);
    }
  }
  let point = new PointG1(new Fp(x), y);
  return point.clearCofactor();
}
