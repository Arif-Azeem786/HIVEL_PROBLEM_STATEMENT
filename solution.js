// solution.js
// Node.js 16+ recommended
const fs = require('fs');

// ---------- BigInt helpers ----------
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

// ---------- parse base-n string -> BigInt ----------
function parseBigIntBase(valueStr, baseNum) {
  const base = BigInt(baseNum);
  let acc = 0n;
  const s = valueStr.trim().toLowerCase();
  for (let ch of s) {
    let digit;
    if (ch >= '0' && ch <= '9') digit = BigInt(ch.charCodeAt(0) - '0'.charCodeAt(0));
    else if (ch >= 'a' && ch <= 'z') digit = BigInt(ch.charCodeAt(0) - 'a'.charCodeAt(0) + 10);
    else {
      // ignore unexpected characters/spaces
      continue;
    }
    if (digit >= base) {
      throw new Error(`Digit ${ch} >= base ${baseNum}`);
    }
    acc = acc * base + digit;
  }
  return acc;
}

// ---------- read input.json ----------
const raw = fs.readFileSync('input.json', 'utf8');
const data = JSON.parse(raw);

// Extract keys info
const n = Number(data.keys.n);
const k = Number(data.keys.k);
if (!n || !k) {
  console.error("Invalid keys.n or keys.k in input.json");
  process.exit(1);
}

// collect points: keys like "1","2",...
const points = [];
for (const key of Object.keys(data)) {
  if (key === 'keys') continue;
  if (!/^\d+$/.test(key)) continue;
  const obj = data[key];
  const x = BigInt(Number(key)); // x coordinate is numeric index
  const base = Number(obj.base);
  const valStr = String(obj.value);
  const y = parseBigIntBase(valStr, base);
  points.push({ x, y, idx: Number(key) });
}

// sort by x ascending and take first k points
points.sort((a,b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));
if (points.length < k) {
  console.error("Not enough points provided to solve for degree (need k points).");
  process.exit(1);
}
const sel = points.slice(0, k);

// ---------- Lagrange interpolation at x=0 (exact BigInt rational) ----------
//
// C = sum_i y_i * prod_{j!=i} (-x_j) / prod_{j!=i} (x_i - x_j)
// We'll compute each term as numerator/denominator BigInts, simplify per-term by gcd to keep numbers smaller,
// then sum into running rational (num/den) using gcd reductions.
let totalNum = 0n;
let totalDen = 1n;

for (let i = 0; i < sel.length; ++i) {
  const yi = sel[i].y;
  let num = yi; // BigInt
  let den = 1n;
  for (let j = 0; j < sel.length; ++j) {
    if (i === j) continue;
    const xj = sel[j].x;
    const xi = sel[i].x;
    num = num * (-xj);
    den = den * (xi - xj);
    // reduce per step to avoid huge blowup
    const g = bigintGCD(num < 0n ? -num : num, den < 0n ? -den : den);
    if (g > 1n) {
      num /= g;
      den /= g;
    }
  }
  // sum totalNum/totalDen + num/den
  const newNum = totalNum * den + num * totalDen;
  const newDen = totalDen * den;
  const g2 = bigintGCD(newNum < 0n ? -newNum : newNum, newDen < 0n ? -newDen : newDen);
  totalNum = newNum / g2;
  totalDen = newDen / g2;
}

// simplify final fraction
const gFinal = bigintGCD(totalNum < 0n ? -totalNum : totalNum, totalDen < 0n ? -totalDen : totalDen);
totalNum /= gFinal;
totalDen /= gFinal;

// Output
if (totalDen === 1n) {
  console.log(totalNum.toString());
} else if (totalNum % totalDen === 0n) {
  console.log((totalNum / totalDen).toString());
} else {
  console.log(totalNum.toString() + "/" + totalDen.toString());
}
