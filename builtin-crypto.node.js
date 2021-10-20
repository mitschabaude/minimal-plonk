import crypto from "node:crypto";

export { randomBytes, sha512 };

function randomBytes(n) {
  return new Uint8Array(crypto.randomBytes(n));
}

// async to match Web Crypto signature
async function sha512(msg) {
  return new Uint8Array(crypto.createHash("sha512").update(msg).digest());
}
