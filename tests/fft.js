import test from "ava";
import {
  evalPoly,
  evalPolyBarycentric,
  evalPolyFFT,
  evalPolyLagrange,
  interpolateIFFT,
  padLength,
  vectorDiv,
  vectorMod,
  vectorMul,
} from "../src/polynomials.js";
import { batchInverse, mod, modInverse } from "../src/modular-arithmetic.js";
import {
  getAllRootsOfUnity,
  largePrimes,
  randomBigIntLength,
  randomRootOfUnity,
} from "../src/random-primes.js";

const p = largePrimes[256];
let randomElement = () => randomBigIntLength(32, false) % p;
let randomVector = (n) => Array(n).fill(0n).map(randomElement);

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

  let f1 = interpolateIFFT(Ff, w, p);
  assert(arrayEqual(f, f1));
});

test("roots-of-unity", ({ assert, is }) => {
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
  let k = 10;
  let n = 2 ** k;
  let W = getAllRootsOfUnity(k, p);
  is(W.length, n, "root is primitive");

  // create polynomial, padd to roots of unity length
  let N = 50;
  let f = Array(N)
    .fill(0)
    .map((_, i) => BigInt(i));

  // fft, ifft should recover original polynomial (up to padded zeros)
  let Ff = evalPolyFFT(f, W, p);
  is(Ff.length, n);

  let f1 = interpolateIFFT(Ff, W, p);
  assert(arrayEqual(f, f1.slice(0, N)));
  assert(arrayEqual(padLength(f, n), f1));

  // explicitly check that f(w^0) === f0 + f1 + ...
  let sumf = f.reduce((sum, fi) => sum + fi, 0n);
  is(sumf, Ff[0]);
});

test("batch-inverse", ({ assert, is }) => {
  let a = randomVector(10);
  let ainv = batchInverse(a, p);
  let ainv1 = a.map((ai) => modInverse(ai, p));
  is(ainv[7], ainv1[7]);
  assert(arrayEqual(ainv, ainv1));
});

test("polynomial-algebra", ({ assert }) => {
  let f = [2n, 1n]; // x + 2
  let g = vectorMod([-2n, 1n], p); // x - 2
  let fg = vectorMod([-4n, 0n, 1n], p); // (x + 2)(x - 2) = x^2 - 4
  let ff = [4n, 4n, 1n]; // (x + 2)(x + 2) = x^2 + 4x + 4

  let k = 5;
  let n = 1 << k;
  let W = getAllRootsOfUnity(k, p);

  let Ff = evalPolyFFT(f, W, p);
  let Fg = evalPolyFFT(g, W, p);
  let Ffg = vectorMul(Ff, Fg, p);
  let Fff = vectorMul(Ff, Ff, p);

  // polynomial multiplication
  let fg0 = interpolateIFFT(Ffg, W, p);
  let ff0 = interpolateIFFT(Fff, W, p);

  assert(arrayEqual(padLength(fg, n), fg0));
  assert(arrayEqual(padLength(ff, n), ff0));

  // polynomial division
  let f0 = interpolateIFFT(vectorDiv(Fff, Ff, p), W, p);
  let g0 = interpolateIFFT(vectorDiv(Ffg, Ff, p), W, p);

  assert(arrayEqual(padLength(f, n), f0));
  assert(arrayEqual(padLength(g, n), g0));

  // map f to (x => f(w*x)) ~= [a0, w*a1, w^2*a2, ...] "shift"
  // is the same as left-shift in evaluation space, because:
  // [f(w*1), f(w*w), ..., f(w*w^(n-1))] = [f(w), f(w^2), ..., f(1)]
  let fshift = vectorMul(f, W, p);
  let fshift0 = interpolateIFFT([...Ff.slice(1), Ff[0]], W, p);
  assert(arrayEqual(fshift, fshift0));
});

test("polynomial-evaluation", ({ is }) => {
  let f = [45n, 0n, 4n, 0n, 1n]; // 45 + 4x^2 + x^4
  let z = 3n;
  let fz = 162n; // f(3) = 5*9 + 4*3^2 + 3^4 = 2 * 9*9

  let fz0 = evalPoly(f, z, p);
  is(fz, fz0);

  let W = getAllRootsOfUnity(5, p);

  let Ff = evalPolyFFT(f, W, p);
  let fz1 = evalPolyBarycentric(Ff, z, W, p);
  let fz2 = evalPolyLagrange(Ff, z, W, p);
  is(fz, fz1);
  is(fz, fz2);
});

function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}
