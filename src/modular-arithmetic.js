// all functions operate on BigInt
export { mod, modExp, modExpNoPrime, modInvert };

// fast modular exponentiation, a^n % p
function modExp(a, n, p) {
  a = mod(a, p);
  // this assumes that p is prime, so that a^(p-1) % p = 1
  n = mod(n, p - 1n);
  let x = 1n;
  for (; n > 0n; n >>= 1n) {
    if (n & 1n) x = mod(x * a, p);
    a = mod(a * a, p);
  }
  return x;
}
// fast modular exponentiation, a^n % q, whithout assuming q is prime
function modExpNoPrime(a, n, q) {
  a = mod(a, q);
  if (n < 0n) throw Error("n must be positive");
  let x = 1n;
  for (; n > 0n; n >>= 1n) {
    if (n & 1n) x = mod(x * a, q);
    a = mod(a * a, q);
  }
  return x;
}

// inverting with EGCD, 1/a in Z_p
function modInvert(a, p) {
  if (a === 0n) throw Error("cannot invert 0");
  a = mod(a, p);
  let b = p;
  let x = 0n;
  let y = 1n;
  let u = 1n;
  let v = 0n;
  while (a !== 0n) {
    let q = b / a;
    let r = mod(b, a);
    let m = x - u * q;
    let n = y - v * q;
    b = a;
    a = r;
    x = u;
    y = v;
    u = m;
    v = n;
  }
  if (b !== 1n) throw Error("inverting failed (no inverse)");
  return mod(x, p);
}

function mod(x, p) {
  x = x % p;
  return x < 0n ? x + p : x;
}