import basis from "../src/basis-256-10.js";
import {
  divideByVanishing,
  evalPolyFFT,
  evalPolyLagrange,
  interpolateIFFT,
  leftShift,
  nextPower2,
  padLength,
  padPermutation,
} from "../src/polynomials.js";
import {
  commit,
  hashTranscript,
  proveEval,
  validateEval,
} from "../src/inner-product.js";
import { mod, modExp, modInverse } from "../src/modular-arithmetic.js";

const FFT_EXPANSION_FACTOR = 4;
const p = basis.p;

export { plonkProve, plonkVerify };

async function plonkProve(circuit, witness) {
  let { n, numberOfColumns, roots, selectors, idPoly, sigmaPoly } =
    plonkSetupCircuit(circuit);
  let transcript = [];

  // round 1 - commit witnesses
  witness = padMany(witness, n);
  let witnessCoeff = witness.map((a) => interpolateIFFT(a, roots, p));
  let Cwitness = witnessCoeff.map((a) => commit(a));
  transcript.push(...Cwitness);

  // round 2 - build permutation check polynomial z(X)
  let beta = await hashTranscript([...transcript, 0n]);
  let gamma = await hashTranscript([...transcript, 1n]);

  let z = Array(n);
  let zCumProd = 1n;
  // TODO batched inverse for O(N) instead of O(N log N)
  for (let i = 0; true; i++) {
    z[i] = zCumProd;
    if (i === n - 1) break;
    let [num, denom] = [1n, 1n];
    for (let j = 0; j < numberOfColumns; j++) {
      num *= witness[j][i] + beta * idPoly[j][i] + gamma;
      denom *= witness[j][i] + beta * sigmaPoly[j][i] + gamma;
    }
    zCumProd = mod(zCumProd * num * modInverse(denom, p), p);
  }
  let zCoeff = interpolateIFFT(z, roots, p);
  let Cz = commit(zCoeff);
  transcript.push(Cz);

  // round 3 - quotient polynomial
  let alpha = await hashTranscript(transcript);
  let factor = FFT_EXPANSION_FACTOR;
  let [expandedRoots] = selectRootsSubset(factor * n, numberOfColumns, basis);

  let fftExpand = (evals) => {
    let coeffs = interpolateIFFT(evals, roots, p);
    let expandedCoeffs = padLength(coeffs, factor * n, 0n);
    return evalPolyFFT(expandedCoeffs, expandedRoots, p);
  };
  let fftExpandMany = (args) => {
    return args.map(fftExpand);
  };

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

  let fullEqCoeff = interpolateIFFT(fullEq, expandedRoots, p);

  let [quotientCoeff, remainder] = divideByVanishing(fullEqCoeff, n, p);
  console.assert(remainder.every((r) => r === 0n));

  let tloCoeff = quotientCoeff.slice(0, n);
  let tmiCoeff = quotientCoeff.slice(n, 2 * n);
  let thiCoeff = quotientCoeff.slice(2 * n, 3 * n);
  let Cquotient = [tloCoeff, tmiCoeff, thiCoeff].map((t) => commit(t));
  transcript.push(...Cquotient);

  // round 4/5 - evaluation
  // TODO tons of efficiency improvements
  let zeta = await hashTranscript(transcript);
  let witnessEval = await Promise.all(
    witnessCoeff.map((a) => proveEval(a, zeta))
  );

  let zEval = await proveEval(zCoeff, zeta);
  let quotientEval = await Promise.all(
    [tloCoeff, tmiCoeff, thiCoeff].map((a) => proveEval(a, zeta))
  );

  let zetaShifted = mod(zeta * roots[1], p);
  let zShiftedEval = await proveEval(zCoeff, zetaShifted);

  let snark = {
    transcript,
    witnessEval,
    zEval,
    zShiftedEval,
    quotientEval,
  };
  return snark;
}

async function plonkVerify(circuit, snark) {
  let { n, numberOfColumns, roots, selectors, sigmaPoly, cofactors } =
    plonkSetupCircuit(circuit);

  // the verifier
  let { transcript, witnessEval, zEval, zShiftedEval, quotientEval } = snark;

  // recompute challenges
  let transcript1 = transcript.slice(0, numberOfColumns);
  let beta = await hashTranscript([...transcript1, 0n]);
  let gamma = await hashTranscript([...transcript1, 1n]);
  let alpha = await hashTranscript(transcript.slice(0, numberOfColumns + 1));
  let zeta = await hashTranscript(transcript);

  // verify evaluations at zeta
  let commits = transcript;
  let evals = [...witnessEval, zEval, ...quotientEval];
  let oks = await Promise.all(
    evals.map((e, i) => validateEval(commits[i], zeta, e.fz, e.proof))
  );
  if (oks.length !== 4 + numberOfColumns || !oks.every((x) => x === true))
    return false;

  // verify evaluation of z at zeta * root
  let zetaShifted = mod(zeta * roots[1], p);
  let ok = await validateEval(
    commits[numberOfColumns],
    zetaShifted,
    zShiftedEval.fz,
    zShiftedEval.proof
  );
  if (ok !== true) return false;

  // evaluate zeta^n - 1, L_0(zeta), selectors, permutations
  let zetaPowerN = modExp(zeta, BigInt(n), p);
  let zH = mod(zetaPowerN - 1n, p);
  let L0 = evalPolyLagrange([1n, ...Array(n - 1).fill(0n)], zeta, roots, p);
  let selectorsAtZeta = selectors.map((q) =>
    evalPolyLagrange(q, zeta, roots, p)
  );
  let sigmaAtZeta = sigmaPoly.map((q) => evalPolyLagrange(q, zeta, roots, p));

  // combine everything to the big equation check
  let [a, b, c] = witnessEval.map((a) => a.fz);
  let z = zEval.fz;
  let zShifted = zShiftedEval.fz;

  let [qlZeta, qrZeta, qoZeta, qmZeta, qcZeta] = selectorsAtZeta;
  let [sigma1, sigma2, sigma3] = sigmaAtZeta;
  let [, k1, k2] = cofactors;
  let [tlo, tmi, thi] = quotientEval.map((t) => t.fz);

  let eq1 = a * b * qmZeta + a * qlZeta + b * qrZeta + c * qoZeta + qcZeta;
  let eq2 =
    z *
      (a + beta * zeta + gamma) *
      (b + beta * k1 * zeta + gamma) *
      (c + beta * k2 * zeta + gamma) -
    zShifted *
      (a + beta * sigma1 + gamma) *
      (b + beta * sigma2 + gamma) *
      (c + beta * sigma3 + gamma);
  let eq3 = (z - 1n) * L0;
  let lhs = mod(eq1 + alpha * eq2 + alpha ** 2n * eq3, p);
  let rhs = mod(zH * (tlo + zetaPowerN * tmi + zetaPowerN ** 2n * thi), p);

  return lhs === rhs;
}

function plonkSetupCircuit({
  selectors,
  permutation,
  circuitLength,
  numberOfColumns,
}) {
  circuitLength ??= selectors[0].length;
  numberOfColumns ??= permutation.length;
  let n = nextPower2(circuitLength);

  let cosets = selectRootsSubset(n, numberOfColumns, basis);
  let cofactors = basis.cofactors.slice(0, numberOfColumns);
  let coroots = cosets.flat();
  let [roots] = cosets;

  selectors = padMany(selectors, n);

  // permutation
  let sigma = padManyPermutations(permutation, n);
  let id = identityPermutation(n, numberOfColumns);

  let permsToPolys = (p) => p.map((s) => s.map((i) => coroots[i]));
  let idPoly = permsToPolys(id);
  let sigmaPoly = permsToPolys(sigma);

  return {
    n,
    numberOfColumns,
    cofactors,
    coroots,
    roots,
    selectors,
    idPoly,
    sigmaPoly,
  };
}

function combinePointwise(combine, ...args) {
  let n = args[0][0].length;
  let combined = Array(n);
  for (let i = 0; i < n; i++) {
    let argsi = args.map((arg) => arg.map((a) => a[i]));
    combined[i] = combine(...argsi, i);
  }
  return combined;
}

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

// function ZH(n, roots, p) {
//   n = BigInt(n);
//   let N = roots.length;
//   let zh = Array(N);
//   for (let i = 0; i < N; i++) {
//     zh[i] = modExp(roots[i], n, p) - 1n;
//   }
//   return zh;
// }
