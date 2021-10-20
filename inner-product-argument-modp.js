// polynomial commitment scheme using the "inner product argument" in a normal prime field Z_p
import basis from "./basis-1024-1e+5.js";
import { sha512 } from "./builtin-crypto.node.js";
import { modExp } from "./modular-arithmetic.js";
const { p, byteLength } = basis;
BigInt.prototype.toJSON = function () {
  return Buffer.from(bigIntArrayToBytes([this], byteLength)).toString("base64");
};

export { commit, proveEval, validateEval };

let f = Array(1000)
  .fill(0)
  .map((_, i) => BigInt(i) ** 2n);
// let f = [-1n, 0n, 1n];
let z = 5;
let comf = commit(f);
let { fz, proof } = await proveEval(f, z);
console.log("proof length:", JSON.stringify([comf, proof, z, fz]).length);
console.log("poly length:", JSON.stringify(f).length);

let ok = await validateEval(comf, z, fz, proof);
console.log("ok", ok);

// console.log(padPower2(f));

function commit(f) {
  return scalarProdGroup(f, basis.G);
}

async function proveEval(f, z) {
  let a = padPower2(f); // ensure length of a is 2^k
  let length = a.length;
  let length0 = length;
  let G = basis.G.slice(0, length);
  let b = zzz(z, length);
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
  let G = basis.G.slice(0, length);
  let b = zzz(z, length);
  let A = comf;
  let ab = fz;
  for (
    let i = 0, half = length >> 1;
    length > 1;
    i++, length = half, half >>= 1
  ) {
    let blo = b.slice(0, half);
    let bhi = b.slice(half, length);
    let Glo = G.slice(0, half);
    let Ghi = G.slice(half, length);
    let [LA, RA, Lab, Rab] = transcript.slice(4 * i, 4 * i + 4);
    let x = await hashTranscript(transcript.slice(0, 4 * i + 4));
    let x2 = multField(x, x);

    b = vecAddField(scalarMultFieldVec(x, blo), bhi);
    G = vecAddGroup(scalarMultGroupVec(x, Glo), Ghi);
    A = (scalarMultGroup(x, A) * LA * scalarMultGroup(x2, RA)) % p;
    ab = (multField(x, ab) + Lab + multField(x2, Rab)) % (p - 1n);
  }
  [G] = G;
  [b] = b;
  return (
    (A - scalarMultGroup(a, G)) % p === 0n &&
    (ab - multField(a, b)) % (p - 1n) === 0n
  );
}

// operations on vectors of scalar field elements / group elements
// both types of vectors are represented as Array<bigint>

// <field[], field[]>
function scalarProdField(f, g, n) {
  if (n === undefined) n = f.length;
  let sum = 0n;
  for (let i = 0; i < n; i++) {
    sum = (sum + f[i] * g[i]) % (p - 1n);
  }
  return sum;
}
// <field[], group[]>
function scalarProdGroup(f, G, n) {
  if (n === undefined) n = f.length;
  let sum = 1n;
  for (let i = 0; i < n; i++) {
    sum = (sum * modExp(G[i], f[i], p)) % p;
  }
  return sum;
}
// field * field
function multField(x, f) {
  return (x * f) % (p - 1n);
}
// field * group
function scalarMultGroup(x, G) {
  return modExp(G, x, p);
}
// field * field[]
function scalarMultFieldVec(x, f) {
  return f.map((fi) => (x * fi) % (p - 1n));
}
// field * group[]
function scalarMultGroupVec(x, G) {
  return G.map((Gi) => modExp(Gi, x, p));
}

// field[] + field[]
function vecAddField(f, g) {
  let n = f.length;
  let h = Array(n);
  for (let i = 0; i < n; i++) {
    h[i] = (f[i] + g[i]) % (p - 1n);
  }
  return h;
}
// group[] + group[]
function vecAddGroup(F, G) {
  let n = F.length;
  let H = Array(n);
  for (let i = 0; i < n; i++) {
    H[i] = (F[i] * G[i]) % p;
  }
  return H;
}

async function hashTranscript(transcript) {
  // sha512 produces 512 bits - is this ok to use as fiat-shamir-random field element?
  return bytesToBigInt(
    await sha512(bigIntArrayToBytes(transcript, byteLength), byteLength)
  );
}

function padPower2(f) {
  // round to next power of 2
  let k = f.length.toString(2).length;
  if (f.length > 1 << k) k++;
  return f.concat(Array((1 << k) - f.length).fill(0n));
}

function zzz(z, n) {
  z = BigInt(z);
  let Z = Array(n);
  let zz = 1n;
  for (let i = 0; i < n; i++) {
    Z[i] = zz;
    zz *= z;
  }
  return Z;
}

function bigIntArrayToBytes(arr, bytesPerBigInt) {
  let n = arr.length;
  let totalBytes = n * bytesPerBigInt;
  let bytes = new Uint8Array(totalBytes);
  for (let i = 0, j = 0; i < n; i++, j += bytesPerBigInt) {
    let bigInt = arr[i];
    for (let k = 0; k < bytesPerBigInt; k++, bigInt >>= 8n) {
      bytes[j + k] = Number(bigInt & 255n);
    }
  }
  return bytes;
}
function bytesToBigInt(bytes) {
  let n = BigInt(bytes.byteLength);
  let bigInt = 0n;
  for (let i = 0n; i < n; i++) {
    bigInt += BigInt(bytes[i]) << (i * 8n);
  }
  return bigInt;
}
