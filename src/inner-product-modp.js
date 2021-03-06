// polynomial commitment scheme using the "inner product argument" in a normal prime field Z_p
import basis from "./basis-modp-256-10.js";
import { bigIntArrayToBytes, bytesToBigInt } from "./bigint.js";
import { sha512 } from "#builtin-crypto";
import { mod, modExp } from "./modular-arithmetic.js";
import { randomBigIntLength } from "./random-primes.js";
import { padPower2 } from "./polynomials.js";
const { p } = basis;

export { commit, proveEval, validateEval, randomFieldElement, hashTranscript };

function commit(f) {
  return scalarProdGroup(f, basis.G);
}

async function proveEval(f, z) {
  let a = padPower2(f); // ensure length of a is 2^k
  let length = a.length;
  let length0 = length;
  let G = basis.G.slice(0, length);
  let b = powersOfZ(z, length);
  let fz = scalarProdField(a, b);
  let transcript = [];

  for (let half = length >> 1; length > 1; length = half, half >>= 1) {
    let alo = a.slice(0, half);
    let ahi = a.slice(half, length);
    let blo = b.slice(0, half);
    let bhi = b.slice(half, length);

    let Glo = G.slice(0, half);
    let Ghi = G.slice(half, length);

    let LA = scalarProdGroup(alo, Ghi);
    let RA = scalarProdGroup(ahi, Glo);
    let Lab = scalarProdField(alo, bhi);
    let Rab = scalarProdField(ahi, blo);

    transcript.push(LA, RA, Lab, Rab);
    let x = await hashTranscript(transcript);

    a = vecAddField(alo, scalarMultFieldVec(x, ahi));
    b = vecAddField(scalarMultFieldVec(x, blo), bhi);
    G = vecAddGroup(scalarMultGroupVec(x, Glo), Ghi);
  }
  [a] = a;
  return { fz, proof: { a, transcript, length: length0 } };
}

async function validateEval(comf, z, fz, proof) {
  let { a, transcript, length } = proof;
  let A = comf;
  let v = fz;
  let xProd = Array(length).fill(1n);
  for (let i = 0, half = length >> 1; half > 0; i++, half >>= 1) {
    let [LA, RA, Lab, Rab] = transcript.slice(4 * i, 4 * i + 4);
    let x = await hashTranscript(transcript.slice(0, 4 * i + 4));
    let x2 = multField(x, x);
    for (let j = 0; j < length; j++) {
      if (!(half & j)) xProd[j] = multField(x, xProd[j]);
    }
    A = addGroup(scalarMultGroup(x, A), LA, scalarMultGroup(x2, RA));
    v = addField(multField(x, v), Lab, multField(x2, Rab));
  }
  let G = scalarProdGroup(xProd, basis.G.slice(0, length));
  let b = scalarProdField(xProd, powersOfZ(z, length));
  let aG = scalarMultGroup(a, G);
  let ab = multField(a, b);
  return A === aG && v === ab;
}

// field/group operations
// elements are BigInts

// <field[], field[]>
function scalarProdField(f, g, n) {
  if (n === undefined) n = f.length;
  let sum = 0n;
  for (let i = 0; i < n; i++) {
    sum = mod(sum + f[i] * g[i], p - 1n);
  }
  return sum;
}
// <field[], group[]>
function scalarProdGroup(f, G, n) {
  if (n === undefined) n = f.length;
  let sum = 1n;
  for (let i = 0; i < n; i++) {
    sum = mod(sum * modExp(G[i], f[i], p), p);
  }
  return sum;
}
// field * field
function multField(x, f) {
  return mod(x * f, p - 1n);
}
// field * group
function scalarMultGroup(x, G) {
  return modExp(G, x, p);
}
// field * field[]
function scalarMultFieldVec(x, f) {
  return f.map((fi) => mod(x * fi, p - 1n));
}
// field * group[]
function scalarMultGroupVec(x, G) {
  return G.map((Gi) => modExp(Gi, x, p));
}
// field + field + ...
function addField(...fs) {
  let sum = 0n;
  for (let f of fs) {
    sum = mod(sum + f, p - 1n);
  }
  return sum;
}
// group + group + ...
function addGroup(...Fs) {
  let sum = 1n;
  for (let F of Fs) {
    sum = mod(sum * F, p);
  }
  return sum;
}
// field[] + field[]
function vecAddField(f, g) {
  let n = f.length;
  let h = Array(n);
  for (let i = 0; i < n; i++) {
    h[i] = mod(f[i] + g[i], p - 1n);
  }
  return h;
}
// group[] + group[]
function vecAddGroup(F, G) {
  let n = F.length;
  let H = Array(n);
  for (let i = 0; i < n; i++) {
    H[i] = mod(F[i] * G[i], p);
  }
  return H;
}
async function hashTranscript(transcript) {
  // sha512 produces 512 bits - is this ok to use as fiat-shamir-random field element?
  let k = basis.byteLength;
  return mod(
    bytesToBigInt(await sha512(bigIntArrayToBytes(transcript, k), k)),
    p - 1n
  );
}
function powersOfZ(z, n) {
  z = BigInt(z);
  let Z = Array(n);
  let zz = 1n;
  for (let i = 0; i < n; i++) {
    Z[i] = mod(zz, p - 1n);
    zz *= z;
  }
  return Z;
}

function randomFieldElement() {
  return randomBigIntLength(basis.byteLength) % p;
}
