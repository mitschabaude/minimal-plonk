export { bigIntArrayToBytes, bytesToBigInt };

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
