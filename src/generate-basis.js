// node generate-basis.js --pbits 256 --dbits 16
import fs from "node:fs";
import minimist from "minimist";
import {
  getAllRootsOfUnity,
  getCosets,
  largePrimes,
  randomBigIntRange,
  randomLargePrime,
} from "./random-primes.js";
let argv = minimist(process.argv.slice(2));

let bitLength = argv.pbits ?? 256;
let degreeBitLength = argv.dbits ?? 10;
let maxColumns = 10;

BigInt.prototype.toJSON = function () {
  return this.toString();
};

let p = largePrimes[bitLength];
if (!p) p = randomLargePrime(bitLength >> 3);
let maxDegree = 1 << degreeBitLength;
console.log(bitLength, maxDegree);

let W = getAllRootsOfUnity(degreeBitLength, p);
let { cosets, cofactors } = getCosets(W, p, maxColumns);

// large collection of random Z_p elements for Pedersen hashing
let G = randomBasisModP(p, maxDegree);

fs.writeFileSync(
  `src/basis-${bitLength}-${degreeBitLength}.js`,
  "let basis = " +
    JSON.stringify({
      bitLength,
      byteLength: bitLength >> 3,
      degreeBitLength,
      maxDegree,
      maxColumns,
      p,
      G,
      W,
      cosets,
      cofactors,
    }) +
    `

let scalars = ["p"];
let arrays = ["G", "W", "cofactors"];
let arrayOfArrays = ["cosets"];
for (let key of scalars) {
  basis[key] = BigInt(basis[key]); 
}
for (let key of arrays) {
  basis[key] = basis[key].map(BigInt); 
}
for (let key of arrayOfArrays) {
  basis[key] = basis[key].map(a => a.map(BigInt)); 
}

export default basis;
`
);

function randomBasisModP(p, length) {
  return Array(length)
    .fill(0n)
    .map(() => randomBigIntRange(0n, p - 1n));
}
