// all functions operate on BigInt
export { mod, modExp, modExpNoPrime, modInverse, batchInverse };

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
function modInverse(a, p) {
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

// a = [a_i] -> [a_i^(-1)] with only (3 + o(1)) multiplications per element
function batchInverse(a, p) {
  let n = a.length;
  // prods = [1, a0, a0*a1, ..., a0*....*a(n-2)], prod = a0*....*a(n-1)
  let prods = Array(n);
  let prod = 1n;
  for (let i = 0; i < n; i++) {
    prods[i] = prod;
    prod = mod(prod * a[i], p);
  }
  // invprod = 1/(a0*....*a(n-1))
  let invprod = modInverse(prod, p);
  // invprods = [1/a0, 1/(a0*a1), ..., 1/(a0*....*a(n-1))]
  let invprods = Array(n);
  for (let i = n - 1; i >= 0; i--) {
    invprods[i] = invprod;
    invprod = mod(a[i] * invprod, p);
  }
  // inv = [1/a0, 1/a1, ..., 1/(a(n-1))]
  // via 1/ai = (a0*...*a(i-1)) / (a0*...*a(i-1)*ai)
  let invs = Array(n);
  for (let i = 0; i < n; i++) {
    invs[i] = mod(prods[i] * invprods[i], p);
  }
  return invs;
}

function mod(x, p) {
  x = x % p;
  return x < 0n ? x + p : x;
}
