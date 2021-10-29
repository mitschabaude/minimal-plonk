import test from "ava";
import {
  evalPolyFFT,
  interpolateInverseFFT,
  padLength,
} from "../src/polynomials.js";
import { mod } from "../src/modular-arithmetic.js";
import {
  getAllRootsOfUnity,
  largePrimes,
  randomRootOfUnity,
} from "../src/random-primes.js";

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

test("roots-of-unity-2^10", ({ assert, is }) => {
  let p = largePrimes[256];
  let k = 10;
  let n = 2 ** k;
  let w = randomRootOfUnity(k, p);

  let W = new Set();
  for (let i = 0, wn = w; i < n; i++, wn = mod(wn * w, p)) {
    W.add(wn);
  }

  assert(W.has(1n), "root is a root");
  is(W.size, n, "root is primitive");
});

test("large-fft-with-new-roots-of-unity", ({ assert, is }) => {
  let p = largePrimes[256];
  let k = 10;
  let n = 2 ** k;
  let W = getAllRootsOfUnity(k, p);
  is(W.length, n, "root is primitive");

  // create polynomial of length 500, padded to roots of unity length
  let N = 500;
  let f0 = Array(N)
    .fill(0)
    .map((_, i) => BigInt(i));
  let f = padLength(f0, n);
  is(f.length, n);

  // fft, ifft should recover original polynomial
  let Ff = evalPolyFFT(f, W, p);
  let f1 = interpolateInverseFFT(Ff, W, p);
  assert(arrayEqual(f, f1));
  assert(arrayEqual(f0, f1.slice(0, N)));

  // explicitly check that f(w^0) === f0 + f1 + ...
  let sumf = f.reduce((sum, fi) => sum + fi, 0n);
  is(sumf, Ff[0]);
});

function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}
