import test from "ava";
import { evalPolyFFT, interpolateInverseFFT } from "../src/fft.js";
import { mod } from "../src/modular-arithmetic.js";

test("small-fft", ({ assert }) => {
  let p = 337n;
  let w1 = 85n;
  let w = [1n];
  for (let wn = w1; wn !== 1n; wn = mod(wn * w1, p)) {
    w.push(wn);
  }

  let f = [3n, 1n, 4n, 1n, 5n, 9n, 2n, 6n];
  let Ff = [31n, 70n, 109n, 74n, 334n, 181n, 232n, 4n];

  let Ff1 = evalPolyFFT(f, w, p);
  assert(arrayEqual(Ff, Ff1));

  let f1 = interpolateInverseFFT(Ff, w, p);
  assert(arrayEqual(f, f1));
});

function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}
