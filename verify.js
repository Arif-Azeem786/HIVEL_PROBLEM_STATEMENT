// verify.js
// Node.js 16+
// Usage: node verify.js
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

function parseBigIntBase(valueStr, baseNum) {
  const base = BigInt(baseNum);
  let acc = 0n;
  const s = valueStr.trim().toLowerCase();
  for (let ch of s) {
    let digit;
    if (ch >= '0' && ch <= '9') digit = BigInt(ch.charCodeAt(0) - '0'.charCodeAt(0));
    else if (ch >= 'a' && ch <= 'z') digit = BigInt(ch.charCodeAt(0) - 'a'.charCodeAt(0) + 10);
    else continue;
    if (digit >= base) throw new Error(`Digit ${ch} >= base ${baseNum}`);
    acc = acc * base + digit;
  }
  return acc;
}

// === read input ===
const raw = fs.readFileSync('input.json', 'utf8');
const data = JSON.parse(raw);
const n = Number(data.keys.n);
const k = Number(data.keys.k);

const points = [];
for (const key of Object.keys(data)) {
  if (key === 'keys') continue;
  if (!/^\d+$/.test(key)) continue;
  const obj = data[key];
  const x = BigInt(Number(key));
  const base = Number(obj.base);
  const y = parseBigIntBase(String(obj.value), base);
  points.push({ x, y, key: Number(key) });
}
points.sort((a,b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
const sel = points.slice(0, k);

// === Ask user to paste computed C here ===
// Replace the value below with your computed result string or number.
// Example: const computedC_str = "-6290016743746469796";
const computedC_str = "-6290016743746469796"; // <<--- paste your computed C here (string)
const computedC = BigInt(computedC_str);

// === Verify polynomial reproduces y_i for every selected x_i ===
// We'll compute P(xp) using Lagrange form (exact rationals) and compare to y_p.

function evaluatePolynomialAt(xp) {
  // compute sum_j yj * prod_{m!=j} (xp - xm) / prod_{m!=j} (xj - xm)
  const m = sel.length;
  const nums = [];
  const dens = [];
  for (let j = 0; j < m; ++j) {
    let num = sel[j].y;
    let den = 1n;
    for (let t = 0; t < m; ++t) {
      if (t === j) continue;
      num = num * (xp - sel[t].x);
      den = den * (sel[j].x - sel[t].x);
      const g = bigintGCD(num < 0n ? -num : num, den < 0n ? -den : den);
      if (g > 1n) { num /= g; den /= g; }
    }
    nums.push(num);
    dens.push(den);
  }
  // combine to common denominator
  let commonDen = 1n;
  for (let d of dens) commonDen *= d;
  let combined = 0n;
  for (let j=0;j<m;++j) {
    const factor = commonDen / dens[j];
    combined += nums[j] * factor;
  }
  // combined / commonDen is P(xp)
  const g = bigintGCD(combined < 0n ? -combined : combined, commonDen < 0n ? -commonDen : commonDen);
  return { num: combined / g, den: commonDen / g };
}

// 1) verify P(0) equals computedC
const P0 = evaluatePolynomialAt(0n);
if (P0.den === 1n) {
  if (P0.num === computedC) console.log("VERIFIED: P(0) equals computed C exactly.");
  else console.error("Mismatch at x=0: P(0) =", P0.num.toString(), "computedC =", computedC.toString());
} else {
  // rational - compare computedC * den ?= num
  if (computedC * P0.den === P0.num) console.log("VERIFIED: P(0) equals computed C (rational match).");
  else console.error("Mismatch at x=0 (rational): P(0) =", P0.num.toString()+"/"+P0.den.toString(), "computedC =", computedC.toString());
}

// 2) verify P(xi) == yi for each selected point
let allMatch = true;
for (let idx=0; idx<sel.length; ++idx) {
  const xp = sel[idx].x;
  const yp = sel[idx].y;
  const Px = evaluatePolynomialAt(xp);
  // compare Px to yp
  if (Px.den === 1n) {
    if (Px.num !== yp) { console.error(`Mismatch at x=${xp.toString()}: P(x)=${Px.num.toString()}, y=${yp.toString()}`); allMatch=false; }
  } else {
    if (yp * Px.den !== Px.num) { console.error(`Mismatch at x=${xp.toString()} (rational): P(x)=${Px.num.toString()}/${Px.den.toString()}, y=${yp.toString()}`); allMatch=false; }
  }
}
if (allMatch) console.log("VERIFIED: polynomial reproduces all selected k points exactly.");
