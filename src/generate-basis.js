// node generate-basis.js --dbits 16
import fs from "node:fs";
import minimist from "minimist";
import { getAllRootsOfUnity, getCosets } from "./random-primes.js";
import bls from "@noble/bls12-381";
import { randomCurvePoint } from "./bls12-381.js";
let argv = minimist(process.argv.slice(2));

let bitLength = 256;
let degreeBitLength = argv.dbits ?? 10;
let maxColumns = 10;

BigInt.prototype.toJSON = function () {
  return this.toString();
};

let p = bls.CURVE.r;
let maxDegree = 1 << degreeBitLength;
console.log(bitLength, maxDegree);

let W = getAllRootsOfUnity(degreeBitLength, p);
let { cosets, cofactors } = getCosets(W, p, maxColumns);

// large collection of random curve elements for Pedersen hashing
let G = randomCurveBasis(maxDegree);

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
let arrays = ["W", "cofactors"];
let arrayOfArrays = ["cosets"];
let arrayOfPoints = ["G"];
for (let key of scalars) {
  basis[key] = BigInt(basis[key]); 
}
for (let key of arrays) {
  basis[key] = basis[key].map(BigInt); 
}
for (let key of arrayOfArrays) {
  basis[key] = basis[key].map(a => a.map(BigInt)); 
}
for (let key of arrayOfPoints) {
  basis[key] = basis[key].map(a => {
    a.x = BigInt(a.x);
    a.y = BigInt(a.y);
    return a;
  }); 
}

export default basis;
`
);

function randomCurveBasis(length) {
  return Array(length)
    .fill(0n)
    .map(() => {
      let point = randomCurvePoint();
      let [x, y] = point.toAffine();
      x = x.value;
      y = y.value;
      return { x, y };
    });
}
