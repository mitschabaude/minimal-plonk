import test from "ava";
import basis from "../src/basis-256-10.js";
import {
  getAllRootsOfUnity,
  getCosets,
  largePrimes,
} from "../src/random-primes.js";
import {
  divideByVanishing,
  evalPolyFFT,
  interpolateIFFT,
  leftShift,
  nextPower2,
  padLength,
  padPermutation,
  vectorMod,
  vectorMul,
} from "../src/polynomials.js";
import { commit, hashTranscript } from "../src/inner-product-modp.js";
import { mod, modExp, modInverse } from "../src/modular-arithmetic.js";

const FFT_EXPANSION_FACTOR = 4;

test("cosets", ({ is }) => {
  const p = largePrimes[256];
  let k = 12;
  let n = 1 << k;
  let W = getAllRootsOfUnity(k, p);

  let {
    cosets: [, k1W, k2W],
  } = getCosets(W, p, 2);

  let S = new Set(W.concat(k1W, k2W));
  is(S.size, 3 * n);
});

test("plonk", async ({ assert, is }) => {
  const p = largePrimes[256];
  let transcript = [];

  let circuitLength = 4;
  let numberOfColumns = 3;
  let n = nextPower2(circuitLength);

  let cosets = selectRootsSubset(n, numberOfColumns, basis);
  let coroots = cosets.flat();
  let [roots] = cosets;

  is(new Set(coroots).size, numberOfColumns * n);

  // // selector polynomials from "PLONK by hand"
  let ql = [0n, 0n, 0n, 1n];
  let qr = [0n, 0n, 0n, 1n];
  let qo = vectorMod([-1n, -1n, -1n, -1n], p);
  let qm = [1n, 1n, 1n, 0n];
  let qc = [0n, 0n, 0n, 0n];
  let selectors = ([ql, qr, qo, qm, qc] = padMany([ql, qr, qo, qm, qc], n));

  // permutation
  let sigma = [
    [n, n + 1, n + 2, 2 * n],
    [0, 1, 2, 2 * n + 1],
    [3, n + 3, 2 * n + 3, 2 * n + 2],
  ];
  sigma = padManyPermutations(sigma, n);
  let id = identityPermutation(n, numberOfColumns);

  let permsToPolys = (p) => p.map((s) => s.map((i) => coroots[i]));
  let idPoly = permsToPolys(id);
  let sigmaPoly = permsToPolys(sigma);
  // let s1 = permToPoly(sigma[0]);
  // let s2 = permToPoly(sigma[1]);
  // let s3 = permToPoly(sigma[2]);

  // round 1 - witness values
  let a = [3n, 4n, 5n, 9n];
  let b = [3n, 4n, 5n, 16n];
  let c = [9n, 16n, 25n, 25n];
  let witness = ([a, b, c] = padMany([a, b, c], n));
  let witnessCoeff = witness.map((a) => interpolateIFFT(a, roots, p));
  let Cwitness = witnessCoeff.map((a) => commit(a));
  transcript.push(...Cwitness);

  // round 2 - build permutation check polynomial z(X)
  let beta = await hashTranscript([...transcript, 0n]);
  let gamma = await hashTranscript([...transcript, 1n]);

  let z = Array(n);
  let zCumProd = 1n;
  // TODO batched inverse for O(N) instead of O(N log N)
  for (let i = 0; i < n; i++) {
    z[i] = zCumProd;
    // adding an early break here would be more efficient, but for now we
    // want the last value for the assertment after the loop
    // if (i === n - 1) break;
    let [num, denom] = [1n, 1n];
    for (let j = 0; j < numberOfColumns; j++) {
      num *= witness[j][i] + beta * idPoly[j][i] + gamma;
      denom *= witness[j][i] + beta * sigmaPoly[j][i] + gamma;
    }
    zCumProd = mod(zCumProd * num * modInverse(denom, p), p);
  }
  is(zCumProd, 1n);
  let zCoeff = interpolateIFFT(z, roots, p);
  let Cz = commit(zCoeff);
  transcript.push(Cz);

  // round 3 - quotient polynomial
  let alpha = await hashTranscript(transcript);
  let factor = FFT_EXPANSION_FACTOR;
  let [expandedRoots] = selectRootsSubset(factor * n, numberOfColumns, basis);

  function fftExpand(evals) {
    let coeffs = interpolateIFFT(evals, roots, p);
    let expandedCoeffs = padLength(coeffs, factor * n, 0n);
    return evalPolyFFT(expandedCoeffs, expandedRoots, p);
  }
  function fftExpandMany(args) {
    return args.map(fftExpand);
  }

  let xWitness = fftExpandMany(witness);
  let xSelectors = fftExpandMany(selectors);
  let xId = fftExpandMany(idPoly);
  let xSigma = fftExpandMany(sigmaPoly);

  let zShifted = leftShift(z);
  let L0 = [1n, ...Array(n - 1).fill(0n)];
  let [xZ, xZShifted, xL0] = fftExpandMany([z, zShifted, L0]);

  let gateEq = combinePointwise(
    ([a, b, c], [ql, qr, qo, qm, qc]) => {
      return mod(a * b * qm + a * ql + b * qr + c * qo + qc, p);
    },
    xWitness,
    xSelectors
  );
  let permutationEq = combinePointwise(
    ([z, zShifted], witness, idPoly, sigmaPoly) => {
      let [num, denom] = [1n, 1n];
      for (let j = 0; j < numberOfColumns; j++) {
        num *= witness[j] + beta * idPoly[j] + gamma;
        denom *= witness[j] + beta * sigmaPoly[j] + gamma;
      }
      return mod(z * num - zShifted * denom, p);
    },
    [xZ, xZShifted],
    xWitness,
    xId,
    xSigma
  );
  let permutationStartEq = combinePointwise(
    ([z, L0]) => {
      return mod((z - 1n) * L0, p);
    },
    [xZ, xL0]
  );

  let fullEq = combinePointwise(
    ([e1, e2, e3]) => mod(e1 + alpha * e2 + alpha ** 2n * e3, p),
    [gateEq, permutationEq, permutationStartEq]
  );

  // check that all equations == 0 on the original roots of unity
  for (let i = 0; i < n * factor; i += factor) {
    is(gateEq[i], 0n);
    is(permutationEq[i], 0n);
    is(permutationStartEq[i], 0n);
    is(fullEq[i], 0n);
  }

  let degree = (f) => {
    let i = [...f].reverse().findIndex((x) => x !== 0n);
    if (i === -1) i = f.length;
    return f.length - i - 1;
  };

  let fullEqCoeff = interpolateIFFT(fullEq, expandedRoots, p);

  let [quotientCoeff, remainder] = divideByVanishing(fullEqCoeff, n, p);
  assert(remainder.every((r) => r === 0n));
  // console.log(fullEqCoeff);

  let quotient = evalPolyFFT(quotientCoeff, expandedRoots, p);
  let zH = ZH(n, expandedRoots, p);
  let zHCoeff = interpolateIFFT(zH, expandedRoots, p);

  is(degree(zHCoeff), n);
  is(degree(quotientCoeff), 3 * n - 4);
  assert(arrayEqual(vectorMul(quotient, zH, p), fullEq));

  let tloCoeff = quotientCoeff.slice(0, n);
  let tmiCoeff = quotientCoeff.slice(n, 2 * n);
  let thiCoeff = quotientCoeff.slice(2 * n, 3 * n);
  let Cquotient = [tloCoeff, tmiCoeff, thiCoeff].map((t) => commit(t));
  transcript.push(...Cquotient);
});

function combinePointwise(combine, ...args) {
  let n = args[0][0].length;
  let combined = Array(n);
  for (let i = 0; i < n; i++) {
    let argsi = args.map((arg) => arg.map((a) => a[i]));
    combined[i] = combine(...argsi, i);
  }
  return combined;
}

function ZH(n, roots, p) {
  n = BigInt(n);
  let N = roots.length;
  let zh = Array(N);
  for (let i = 0; i < N; i++) {
    zh[i] = modExp(roots[i], n, p) - 1n;
  }
  return zh;
}

function polyMUltiply(f, g, p) {}

function selectRootsSubset(degree, columns, { cosets, maxDegree, maxColumns }) {
  // given n, finds next power of 2 and selects roots of unity and cosets with that length
  if (nextPower2(degree) !== degree) {
    throw Error("degree must be power of 2");
  }
  if (columns > maxColumns) {
    throw Error(`only ${maxColumns} columns supported, got ${columns}`);
  }
  if (degree > maxDegree) {
    throw Error(`only degree up to ${maxDegree} supported, got ${degree}`);
  }
  let m = maxDegree / degree;
  return cosets.slice(0, columns).map((W) => selectEveryMth(W, m));
}

function selectEveryMth(a, m) {
  return a.filter((_, i) => i % m === 0);
}

function identityPermutation(n, k) {
  return Array(k)
    .fill(0)
    .map((_, j) =>
      Array(n)
        .fill(0)
        .map((_, i) => j * n + i)
    );
}

function padMany(args, n, fillValue = 0n) {
  return args.map((f) => padLength(f, n, fillValue));
}
function padManyPermutations(args, n) {
  return args.map((f, i) => padPermutation(f, n, i * n));
}

function arrayEqual(a, b) {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}
