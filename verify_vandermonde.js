// verify_vandermonde.js
// Node.js 16+
// Exact rational arithmetic using BigInt

const fs = require('fs');

function bigintGCD(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

function normalizeRat(r) {
  // r is {n: BigInt, d: BigInt}
  if (r.d < 0n) { r.n = -r.n; r.d = -r.d; }
  if (r.n === 0n) { r.d = 1n; return r; }
  const g = bigintGCD(r.n < 0n ? -r.n : r.n, r.d);
  r.n /= g; r.d /= g;
  return r;
}

function addRat(a, b) {
  const n = a.n * b.d + b.n * a.d;
  const d = a.d * b.d;
  return normalizeRat({n, d});
}
function subRat(a, b) {
  const n = a.n * b.d - b.n * a.d;
  const d = a.d * b.d;
  return normalizeRat({n, d});
}
function mulRat(a, b) {
  const n = a.n * b.n;
  const d = a.d * b.d;
  return normalizeRat({n, d});
}
function divRat(a, b) {
  if (b.n === 0n) throw new Error("division by zero");
  const n = a.n * b.d;
  const d = a.d * b.n;
  return normalizeRat({n, d});
}
function fromBigInt(x) { return normalizeRat({n: BigInt(x), d: 1n}); }

// parse base-N to BigInt
function parseBigIntBase(valueStr, baseNum) {
  const base = BigInt(baseNum);
  let acc = 0n;
  const s = String(valueStr).trim().toLowerCase();
  for (let ch of s) {
    let digit;
    if (ch >= '0' && ch <= '9') digit = BigInt(ch.charCodeAt(0) - '0'.charCodeAt(0));
    else if (ch >= 'a' && ch <= 'z') digit = BigInt(ch.charCodeAt(0) - 'a'.charCodeAt(0) + 10);
    else { throw new Error("Invalid char in value: "+ch); }
    if (digit >= base) throw new Error("Digit >= base");
    acc = acc * base + digit;
  }
  return acc;
}

// Read input
const data = JSON.parse(fs.readFileSync('input.json','utf8'));
const numeric_keys = Object.keys(data).filter(k=>k!=='keys').map(k=>Number(k)).filter(x=>!isNaN(x)).sort((a,b)=>a-b);
const k = Number(data.keys.k);
if (numeric_keys.length < k) throw new Error("Not enough points");

// select first k
const sel = numeric_keys.slice(0,k);
console.log("Selected keys:", sel.join(", "));

// prepare xs and ys (BigInt)
const xs = sel.map(x => BigInt(x));
const ys = sel.map(key => {
  const obj = data[String(key)];
  return parseBigIntBase(String(obj.value), Number(obj.base));
});

// Build Vandermonde matrix V_ij = x_i^(j) for i=0..k-1, j=0..k-1 (columns: x^0, x^1, ...)
const V = Array.from({length:k}, () => Array(k).fill(null));
for (let i=0;i<k;++i) {
  let power = 1n;
  for (let j=0;j<k;++j) {
    V[i][j] = fromBigInt(power); // rational rep
    power = power * xs[i];
  }
}

// RHS vector of rationals
const b = ys.map(y => fromBigInt(y));

// Gaussian elimination (in-place) on rationals: transform V to upper triangular and solve
function gaussianSolve(A, bvec) {
  const n = A.length;
  // make copies
  const M = A.map(row => row.map(r => ({n:r.n,d:r.d})));
  const B = bvec.map(r => ({n:r.n,d:r.d}));
  for (let col=0; col<n; ++col) {
    // pivot: find row >= col with non-zero in M[row][col]
    let pivot = -1;
    for (let r=col;r<n;++r) {
      if (M[r][col].n !== 0n) { pivot = r; break; }
    }
    if (pivot === -1) throw new Error("Singular matrix (no pivot)");
    // swap rows if needed
    if (pivot !== col) {
      [M[pivot], M[col]] = [M[col], M[pivot]];
      [B[pivot], B[col]] = [B[col], B[pivot]];
    }
    // normalize pivot row: divide entire row by pivot element
    const pivVal = M[col][col];
    for (let c=col;c<n;++c) M[col][c] = divRat(M[col][c], pivVal);
    B[col] = divRat(B[col], pivVal);
    // eliminate below
    for (let r=col+1;r<n;++r) {
      const factor = M[r][col];
      if (factor.n === 0n) continue;
      for (let c=col;c<n;++c) {
        M[r][c] = subRat(M[r][c], mulRat(factor, M[col][c]));
      }
      B[r] = subRat(B[r], mulRat(factor, B[col]));
    }
  }
  // back substitution
  const x = Array(n).fill(fromBigInt(0n));
  for (let i=n-1;i>=0;--i) {
    let sum = fromBigInt(0n);
    for (let j=i+1;j<n;++j) sum = addRat(sum, mulRat(M[i][j], x[j]));
    x[i] = subRat(B[i], sum);
  }
  return x; // array of rationals {n,d}
}

// Solve
const coeffs = gaussianSolve(V,b); // coeffs[0] is a0 (constant)
console.log("Coefficient a0 (rational):", coeffs[0].n.toString()+"/"+coeffs[0].d.toString());
if (coeffs[0].d === 1n) console.log("a0 integer:", coeffs[0].n.toString());

// Convert coeffs to BigInt integer if denominators are 1 (or extract reduced numerator/den)
const allInt = coeffs.every(r => r.d === 1n);
if (!allInt) console.log("Some coefficients are rational, but a0 above is exact.");

// Now compute P(0) via coefficients and compare with Lagrange method
const a0 = coeffs[0];
function evaluatePolyWithCoeffsAt(t) {
  // compute sum_{j} a_j * t^j
  let tPow = 1n;
  let acc = {n:0n,d:1n};
  for (let j=0;j<coeffs.length;++j) {
    const term = mulRat(coeffs[j], fromBigInt(tPow));
    acc = addRat(acc, term);
    tPow *= BigInt(t);
  }
  return normalizeRat(acc);
}

// Lagrange evaluate at t (rational)
function evaluateLagrangeAt(t) {
  // sum_i y_i * product_{j!=i} (t - x_j)/(x_i - x_j)
  let num = {n:0n,d:1n};
  for (let i=0;i<k;++i) {
    let term = fromBigInt(ys[i]);
    for (let j=0;j<k;++j) {
      if (i===j) continue;
      term = mulRat(term, normalizeRat({n:(BigInt(t)-xs[j]), d: (xs[i]-xs[j])}));
    }
    num = addRat(num, term);
  }
  return normalizeRat(num);
}

// Compare at t=0 and t=100 (random)
const checkTs = [0, 100];
for (const t of checkTs) {
  const c1 = evaluatePolyWithCoeffsAt(t);
  const c2 = evaluateLagrangeAt(t);
  console.log(`P(${t}) via coeffs = ${c1.n.toString()}/${c1.d.toString()}; via Lagrange = ${c2.n.toString()}/${c2.d.toString()}`);
  if (c1.n * c2.d === c2.n * c1.d) console.log(`Match at t=${t}`);
  else console.log(`Mismatch at t=${t}`);
}
