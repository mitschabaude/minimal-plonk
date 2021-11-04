// fft / ifft in prime field for moving between
// 1) coefficient representation f(x) = f0 + x*f1 + ... + f(n-1)x^(n-1)
// 2) evaluation representation f(1) = Ff0, f(w) = Ff1, ..., f(w^(n-1)) = Ff(n-1)
import { batchInverse, mod, modExp, modInverse } from "./modular-arithmetic.js";

export {
  vectorMod,
  vectorAdd,
  vectorSub,
  vectorMul,
  vectorDiv,
  leftShift,
  divideByVanishing,
  evalPolyFFT,
  interpolateIFFT,
  padPower2,
  padLength,
  padPermutation,
  evalPoly,
  evalPolyBarycentric,
  evalPolyLagrange,
  nextPower2,
};

// vector operations in evaluation space
// these can be turned into operations in coefficient space by wrapping in FFT, IFFT

function vectorMod(a, p) {
  return a.map((ai) => mod(ai, p));
}

function vectorAdd(a, b, p) {
  let n = Math.max(a.length, b.length);
  let plus = Array(n);
  for (let i = 0; i < n; i++) {
    plus[i] = mod((a[i] ?? 0n) + (b[i] ?? 0n), p);
  }
  return plus;
}
function vectorSub(a, b, p) {
  let n = Math.max(a.length, b.length);
  let sub = Array(n);
  for (let i = 0; i < n; i++) {
    sub[i] = mod((a[i] ?? 0n) - (b[i] ?? 0n), p);
  }
  return sub;
}
function vectorMul(a, b, p) {
  let n = Math.max(a.length, b.length);
  let mul = Array(n);
  for (let i = 0; i < n; i++) {
    mul[i] = mod((a[i] ?? 0n) * (b[i] ?? 0n), p);
  }
  return mul;
}
function vectorDiv(a, b, p) {
  return vectorMul(a, batchInverse(b, p), p);
}

function leftShift(a, k = 1) {
  return [...a.slice(k), ...a.slice(0, k)];
}

// coefficient space

function divideByVanishing(f, n, p) {
  // polynomial division f(X) / (X^n - 1) with remainder
  // very cheap, 0 multiplications
  // strategy:
  // start with q(X) = 0, r(X) = f(X)
  // then start changing q, r while preserving the identity:
  // f(X) = q(X) * (X^n - 1) + r(X)
  // in every step, move highest-degree term of r into the product
  // => r eventually has degree < n and we're done
  let q = Array(f.length).fill(0n);
  let r = [...f];
  for (let i = f.length - 1; i >= n; i--) {
    let leadingCoeff = r[i];
    if (leadingCoeff === 0n) continue;
    r[i] = 0n;
    r[i - n] = mod(r[i - n] + leadingCoeff, p);
    q[i - n] = mod(q[i - n] + leadingCoeff, p);
  }
  return [q, r];
}

// given taylor coefficients of f, returns evaluations on mult. subgroup
function evalPolyFFT(f, w, p) {
  let n = w.length;
  if (f.length > n)
    throw Error("degree must be smaller than number of roots of unity");
  if (f.length < n) f = padLength(f, n, 0n);
  let nHalf = n >> 1;
  if (n === 1) return f; // base case: degree-0 poly

  let aEven = f.filter(filterEven);
  let aOdd = f.filter(filterOdd);
  let wHalf = w.filter(filterEven);
  let fEven = evalPolyFFT(aEven, wHalf, p);
  let fOdd = evalPolyFFT(aOdd, wHalf, p);

  let Ff = Array(n);
  for (let i = 0; i < nHalf; i++) {
    let oddPart = w[i] * fOdd[i];
    let evenPart = fEven[i];
    Ff[i] = mod(evenPart + oddPart, p);
    Ff[i + nHalf] = mod(evenPart - oddPart, p);
  }
  return Ff;
}

// given evaluations on mult. subgroup of f, returns taylor coefficients
function interpolateIFFT(Ff, w, p) {
  let f0 = evalPolyFFT(Ff, w, p);
  let n = w.length;
  let ninv = modInverse(BigInt(n), p);
  let f = [f0[0]].concat(f0.slice(1).reverse());
  for (let i = 0; i < n; i++) {
    f[i] = mod(ninv * f[i], p);
  }
  return f;
}

// evaluate a polynomial at some z from coefficient representation in 2N
function evalPoly(f, z, p) {
  let n = f.length;
  let zpow = 1n;
  let sum = 0n;
  for (let i = 0; i < n; i++) {
    sum = mod(sum + f[i] * zpow, p);
    zpow = mod(zpow * z, p);
  }
  return sum;
}

// evaluate a polynomial at some z from evaluation representation in 5N + o(N)
// doing iFFT + normal evaluation would be O(N log(N)) + 2N per evaluation
function evalPolyBarycentric(Ff, z, W, p) {
  let n = Ff.length;
  console.assert(n === W.length);
  let Wmz = W.map((wi) => wi - z);
  // 3N + o(N)
  let denoms = batchInverse(Wmz, p);
  // 2N
  let sum = 0n;
  for (let i = 0; i < n; i++) {
    sum = mod(sum + Ff[i] * W[i] * denoms[i], p);
  }
  // o(N)
  let nn = BigInt(n);
  let zn = modExp(z, nn, p);
  let ninv = modInverse(nn, p);
  return mod((1n - zn) * ninv * sum, p);
}

// more efficient and simpler lagrange evaluation,
// without using the barycentric formula
// 5N + o(N) = 5N - 3 multiplications + 1 inversion of N
function evalPolyLagrange(Ff, z, W, p) {
  let n = Ff.length;
  // forwardProds = [1, (z - w[0]), ..., (z - w[0])...(z - w[n-2))]
  // backwardProds = [(z - w[1])...(z - w[n-1]), ..., (z - w[n-1]), 1]
  // L_i(z) = (1/n) * w[i] * forwardProds[i] * backwardProds[i]
  let forwardProds = Array(n);
  let backwardProds = Array(n);
  let forwardProd = 1n;
  let backwardProd = 1n;
  // 2N - 4
  for (let i = 0; true; i++) {
    forwardProds[i] = forwardProd;
    backwardProds[n - 1 - i] = backwardProd;
    if (i === n - 1) break;
    forwardProd = mod(forwardProd * (z - W[i]), p);
    backwardProd = mod(backwardProd * (z - W[n - 1 - i]), p);
  }
  // N
  let lagrangeEvals = Array(n);
  for (let i = 0; i < n; i++) {
    lagrangeEvals[i] = mod(forwardProds[i] * backwardProds[i], p);
  }
  // 2N
  let sum = 0n;
  for (let i = 0; i < n; i++) {
    sum = mod(sum + Ff[i] * W[i] * lagrangeEvals[i], p);
  }
  // 1 + inversion
  let nn = BigInt(n);
  let ninv = modInverse(nn, p);
  return mod(ninv * sum, p);
}

function filterEven(_, i) {
  return i % 2 === 0;
}
function filterOdd(_, i) {
  return i % 2 === 1;
}

function nextPower2(n) {
  return 1 << Math.ceil((n - 1).toString(2).length);
}

// pad array with zeros up to the next power of 2
function padPower2(f) {
  let n = nextPower2(f.length);
  return f.concat(Array(n - f.length).fill(0n));
}
// pad array with zeros up to length n
function padLength(f, n, fillValue = 0n) {
  if (n < f.length) throw Error("cannot pad up to length n");
  return f.concat(Array(n - f.length).fill(fillValue));
}
// pad array describing a permutation with identity (=dummy) values
function padPermutation(f, n, valueToAdd = 0) {
  if (n < f.length) throw Error("cannot pad up to length n");
  let padded = [...f];
  for (let i = f.length; i < n; i++) {
    padded[i] = i + valueToAdd;
  }
  return padded;
}
