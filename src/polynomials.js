// fft / ifft in prime field for mobing between
// 1) coefficient representation f(x) = f0 + x*f1 + ... + f(n-1)x^(n-1)
// 2) evaluation representation f(1) = Ff0, f(w) = Ff1, ..., f(w^(n-1)) = Ff(n-1)
import { mod, modInvert } from "./modular-arithmetic.js";

export { evalPolyFFT, interpolateInverseFFT, padPower2, padLength };

// given taylor coefficients of f, returns evaluations on mult. subgroup
function evalPolyFFT(f, w, p) {
  let n = f.length;
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
function interpolateInverseFFT(Ff, w, p) {
  let f0 = evalPolyFFT(Ff, w, p);
  let n = Ff.length;
  let ninv = modInvert(BigInt(n), p);
  let f = [f0[0]].concat(f0.slice(1).reverse());
  for (let i = 0; i < n; i++) {
    f[i] = mod(ninv * f[i], p);
  }
  return f;
}

function filterEven(_, i) {
  return i % 2 === 0;
}
function filterOdd(_, i) {
  return i % 2 === 1;
}

// pad array with zeros up to the next power of 2
function padPower2(f) {
  let k = f.length.toString(2).length;
  if (f.length > 1 << k) k++;
  return f.concat(Array((1 << k) - f.length).fill(0n));
}
// pad array with zeros up to length n
function padLength(f, n, fillValue = 0n) {
  if (n < f.length) throw Error("cannot pad up to length n");
  return f.concat(Array(n - f.length).fill(fillValue));
}
