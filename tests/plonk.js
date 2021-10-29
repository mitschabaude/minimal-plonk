import test from "ava";
import basis from "../src/basis-256-10.js";
import {
  getAllRootsOfUnity,
  getTwoCosets,
  largePrimes,
} from "../src/random-primes.js";
import { padLength, padPermutation, vectorMod } from "../src/polynomials.js";

test("cosets", ({ is }) => {
  const p = largePrimes[256];
  let k = 12;
  let n = 1 << k;
  let W = getAllRootsOfUnity(k, p);

  let [k1W, k2W] = getTwoCosets(W, p);

  let S = new Set(W.concat(k1W, k2W));
  is(S.size, 3 * n);
});

test("prepare-plonk", ({ is }) => {
  const p = largePrimes[256];

  let { W, k1W, k2W, degree: n } = selectRootsSubset(4, basis);
  let Wbig = W.concat(k1W, k2W);

  is(W.length, n);
  is(Wbig.length, 3 * n);

  // // selector polynomials from "PLONK by hand"
  let ql = [0n, 0n, 0n, 1n];
  let qr = [0n, 0n, 0n, 1n];
  let qo = vectorMod([-1n, -1n, -1n, -1n], p);
  let qm = [1n, 1n, 1n, 0n];
  let qc = [0n, 0n, 0n, 0n];

  // permutation
  let N = BigInt(n);
  let sigma = [
    [N, N + 1n, N + 2n, 2n * N],
    [0n, 1n, 2n, 2n * N + 1n],
    [3n, N + 3n, 2n * N + 3n, 2n * N + 2n],
  ];
  let s1 = padPermutation(sigma[0], n).map((i) => Wbig[i]);
  let s2 = padPermutation(sigma[1], n).map((i) => Wbig[i]);
  let s3 = padPermutation(sigma[2], n).map((i) => Wbig[i]);

  // witness values
  let a = [3n, 4n, 5n, 9n];
  let b = [3n, 4n, 5n, 16n];
  let c = [9n, 16n, 25n, 25n];
});

function selectRootsSubset(length, { W, k1W, k2W, maxDegree }) {
  // given n, finds next power of 2 and selects roots of unity and cosets with that length
  let k = Math.ceil((length - 1).toString(2).length);
  let degree = 1 << k;
  let m = maxDegree / degree;
  let selectEveryMth = (_, i) => i % m === 0;
  return {
    W: W.filter(selectEveryMth),
    k1W: k1W.filter(selectEveryMth),
    k2W: k2W.filter(selectEveryMth),
    degree,
  };
}
