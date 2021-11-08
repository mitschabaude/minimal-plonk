export { bigIntArrayToBytes, bytesToBigInt, getByteLength };

function bigIntArrayToBytes(arr, bytesPerBigInt) {
  let n = arr.length;
  let totalBytes = n * bytesPerBigInt;
  let bytes = new Uint8Array(totalBytes);
  for (let i = 0, j = 0; i < n; i++, j += bytesPerBigInt) {
    let bigInt = arr[i];
    if (typeof bigInt !== "bigint") {
      console.log(bigInt);
      throw Error("is not bigint");
    }
    for (let k = 0; k < bytesPerBigInt; k++, bigInt >>= 8n) {
      // if (typeof bigInt !== "bigint") {
      //   throw Error("is not bigint");
      // }
      let number = bigInt & 255n;
      bytes[j + k] = Number(number);
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

// how many bytes a bigint has
function getByteLength(x) {
  return Math.ceil(x.toString(16).length / 2);
}
