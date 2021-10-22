// generate random large primes for using F_p
// import { randomBytes } from "./builtin-crypto.js";
import { randomBytes } from "#builtin-crypto";
import first1000Primes from "./first-1000-primes.js";
import { modExp } from "./modular-arithmetic.js";
const knownPrimes = first1000Primes.map((x) => BigInt(x));

export { largePrimes, randomLargePrime, randomBigIntLength, randomBigIntRange };

// random large primes found with miller-rabin test
const p1024 =
  157504965780504614023334664722541728560981130710635270073862108646389561302305908561100988570988580155311773778914144313059319305341352463134445443240472104178221313037238198668989720238174409922800328630409834284445540338696349835638972137809943783541174306819552187964653302389983399933925385355609237488681n;
const p2048 =
  21916817983981762168038706694528552914931674826094186134115819192797625299650082702951626651467493159408604454270438749729738426067903167309453840325929214072071153251368222218363604931492194088379481066598902964207010534185172561356424950573239390543915154285859794070452984251454174107546333476003605657582574039964236737770073112244876504445617403351815091234589431031366301636549295407000671428717790253439342475993260404024202141305526201910799619023814070666620372565862826345680556243739927570772322783335747732144238832563287108793426842525048450376933782008571774571068346124035131399873030449967779641726477n;
const p3096 =
  76896909810841711131873589708863212466728514712262180123451809243111705665915229552881073032199695994413345502284321061728690203813194584461522862116208164086939660772175859283476397220933733151225115181643572704141586505774303442696388262877795175066038225922119223777324562598157200986681175275639007288173114265469772534668207661198473425026791474958998418597140065271261969194013739708987470715523051579165086169770051892210896180092551191708016082654654263439802435188372255210214786541958924945789072210828430275554147141794149680488763747633242517244766372271198141226030981562407288646815652794904891320941558320075619221906427807405791806180914063380513706774487362047257716437999215648783240374876773872482103693060533090875732469267352734429798372866449396556109926124842939599623849398707290038827899220023614308871978985711916625841856518158900557782393445773530420528689679398991771371819217461980026583341943124331397n;
const largePrimes = {
  [1024]: p1024,
  [2048]: p2048,
  [3096]: p3096,
};

// console.log(randomLargePrime(2048 >> 3));

// returns a random, k-byte prime as BigInt type
function randomLargePrime(byteLength) {
  while (true) {
    let p = randomBigIntLength(byteLength);
    if (millerRabinIsOddPrime(p)) return p;
  }
}

// from wikipedia:

// Input #1: n > 3, an odd integer to be tested for primality
// Input #2: k, the number of rounds of testing to perform
// Output: “composite” if n is found to be composite, “probably prime” otherwise

// write n as 2**r·d + 1 with d odd (by factoring out powers of 2 from n − 1)
// WitnessLoop: repeat k times:
//     pick a random integer a in the range [2, n − 2]
//     x ← a**d mod n
//     if x = 1 or x = n − 1 then
//         continue WitnessLoop
//     repeat r − 1 times:
//         x ← x**2 mod n
//         if x = n − 1 then
//             continue WitnessLoop
//     return “composite”
// return “probably prime”

function millerRabinIsOddPrime(n) {
  const k = 10;
  if (n === 2n || n === 3n) return true;
  if (n < 2n) return false;
  // check if divisible by one of first 1000 primes
  for (let p of knownPrimes) {
    if (n % p === 0n && n > p) return false;
  }

  // write n - 1 = 2^r * d, d odd
  let d = n - 1n;
  let r = 0n;
  for (; d % 2n !== 0n; d /= 2n, r++);

  WitnessLoop: for (let i = 0; i < k; i++) {
    let a = randomBigIntRange(2n, n - 2n);
    // let x = a ** d % n; // too large
    let x = modExp(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    for (let i = 0; i + 1 < r; i++) {
      x = (x * x) % n;
      if (x === n - 1n) continue WitnessLoop;
    }
    return false;
  }
  return true;
}

function randomBigIntRange(min, max) {
  while (true) {
    let n = 0n;
    let length = Math.ceil((max - min).toString(16).length / 2);
    let lengthn = BigInt(length);
    let bytes = randomBytes(length);
    for (let i = 0n; i < lengthn; i++) {
      n += BigInt(bytes[i]) << (8n * i);
    }
    if (n <= max - min) return min + n;
  }
}

function randomBigIntLength(byteLength) {
  let lengthn = BigInt(byteLength);
  let bytes = randomBytes(byteLength);
  let n = 0n;
  for (let i = 0n; i < lengthn; i++) {
    let byte = bytes[i];
    if (i === lengthn - 1n && byte < 128) {
      byte += 128;
    }
    n += BigInt(byte) << (8n * i);
  }
  return n;
}
