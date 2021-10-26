export { randomBytes, sha512, sha256 };

function randomBytes(n) {
  return crypto.getRandomValues(new Uint8Array(n));
}

async function sha512(msg) {
  return new Uint8Array(await crypto.subtle.digest("SHA-512", msg));
}
async function sha256(msg) {
  return new Uint8Array(await crypto.subtle.digest("SHA-256", msg));
}
