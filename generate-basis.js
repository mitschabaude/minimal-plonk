// node generate-basis.js --nbits 1024 --nbasis 10000
import fs from "node:fs";
import minimist from "minimist";
import {
  largePrimes,
  randomBigIntRange,
  randomLargePrime,
} from "./random-primes.js";
let argv = minimist(process.argv.slice(2));

BigInt.prototype.toJSON = function () {
  return this.toString();
};

let bitLength = argv.nbits ?? 1024;
let p = largePrimes[bitLength];
if (!p) p = randomLargePrime(bitLength >> 3);
let n = argv.nbasis ?? 10000;
console.log(bitLength, n);

fs.writeFileSync(
  `basis-${bitLength}-${n.toExponential()}.js`,
  "let basis = " +
    JSON.stringify({
      p,
      bitLength,
      byteLength: bitLength >> 3,
      // large collection of random Z_p elements, e.g. for Pedersen hashing
      G: randomBasisModP(p, n),
      // some additional random z_p elements
      Q: randomBigIntRange(0n, p - 1n),
      R: randomBigIntRange(0n, p - 1n),
    }) +
    `

basis.p = BigInt(basis.p);
basis.Q = BigInt(basis.Q);
basis.R = BigInt(basis.R);
basis.G = basis.G.map((x) => BigInt(x));

export default basis;
`
);

function randomBasisModP(p, length) {
  return Array(length)
    .fill(0n)
    .map(() => randomBigIntRange(0n, p - 1n));
}
