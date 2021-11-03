import test from "ava";
import basis from "../src/basis-256-10.js";
import {
  getAllRootsOfUnity,
  getTwoCosets,
  largePrimes,
} from "../src/random-primes.js";
import { padPermutation, vectorMod } from "../src/polynomials.js";
import { commit, hashTranscript } from "../src/inner-product-modp.js";
import { mod, modInverse } from "../src/modular-arithmetic.js";

test("cosets", ({ is }) => {
  const p = largePrimes[256];
  let k = 12;
  let n = 1 << k;
  let W = getAllRootsOfUnity(k, p);

  let [k1W, k2W] = getTwoCosets(W, p);

  let S = new Set(W.concat(k1W, k2W));
  is(S.size, 3 * n);
});

test("prepare-plonk", async ({ is }) => {
  const p = largePrimes[256];
  let transcript = [];

  let { W, k1W, k2W, degree: n } = selectRootsSubset(4, basis);
  let WBig = W.concat(k1W, k2W);

  is(W.length, n);
  is(WBig.length, 3 * n);

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
  let permToPoly = (s, N) => padPermutation(s, n, N).map((i) => WBig[i]);
  let s1 = permToPoly(sigma[0], 0n);
  let s2 = permToPoly(sigma[1], N);
  let s3 = permToPoly(sigma[2], 2n * N);

  // round 1 - witness values
  let a = [3n, 4n, 5n, 9n];
  let b = [3n, 4n, 5n, 16n];
  let c = [9n, 16n, 25n, 25n];
  let Ca = commit(a);
  let Cb = commit(b);
  let Cc = commit(c);
  transcript.push(Ca, Cb, Cc);

  // round 2 - build permutation check polynomial z(X)
  let beta = await hashTranscript([...transcript, 0n]);
  let gamma = await hashTranscript([...transcript, 1n]);

  let z = Array(n);
  let zCumProd = 1n;
  // TODO batched inverse for O(N) instead of O(N log N)
  for (let i = 0; true; i++) {
    z[i] = zCumProd;
    if (i === n - 1) break;
    let num =
      (a[i] + beta * W[i] + gamma) *
      (b[i] + beta * k1W[i] + gamma) *
      (c[i] + beta * k2W[i] + gamma);
    let denom =
      (a[i] + beta * s1[i] + gamma) *
      (b[i] + beta * s2[i] + gamma) *
      (c[i] + beta * s3[i] + gamma);
    zCumProd = mod(zCumProd * num * modInverse(denom, p), p);
  }
  let Cz = commit(z);
  transcript.push(Cz);

  // round 3 - quotient polynomial
  let alpha = await hashTranscript(transcript);

  let equation1;
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
