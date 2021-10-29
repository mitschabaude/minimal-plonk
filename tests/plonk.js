import test from "ava";
import {
  getAllRootsOfUnity,
  getTwoCosets,
  largePrimes,
} from "../src/random-primes.js";

test("cosets", ({ is }) => {
  const p = largePrimes[256];
  let k = 12;
  let n = 1 << k;
  let W = getAllRootsOfUnity(k, p);

  let [k1W, k2W] = getTwoCosets(W, p);

  let S = new Set(W.concat(k1W, k2W));
  is(S.size, 3 * n);
});
