/**
 * Real Cross-Language Node.js vs Python Parity Test (`tests/parity_js_check.mjs`).
 * Loads `results/model_artifact.json` and `results/parity_150_profiles.json` exported by Python,
 * scores all 150 profiles using the JavaScript engine (`web/js/engine.mjs`), and asserts
 * maximum discrepancy < 1e-7 across all outputs.
 */
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { scoreProfileFromArtifactJS } from '../web/js/engine.mjs';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, '..');

const modelArt = JSON.parse(fs.readFileSync(path.join(repoRoot, 'results', 'model_artifact.json'), 'utf8'));
const vectors = JSON.parse(fs.readFileSync(path.join(repoRoot, 'results', 'parity_150_profiles.json'), 'utf8'));

if (modelArt.n_analyzed !== 4954) {
  throw new Error(`Expected n_analyzed = 4954, got ${modelArt.n_analyzed}`);
}
if (vectors.length !== 150) {
  throw new Error(`Expected 150 parity vectors, got ${vectors.length}`);
}

let maxDiff = 0.0;
for (const vec of vectors) {
  const res = scoreProfileFromArtifactJS(modelArt, vec.profile, 1.0);
  const c36 = res.curve[4];
  const c60 = res.curve[6];
  const c120 = res.curve[9];

  // Partition of unity check
  const unityErr = Math.abs((c60.os + c60.cif1 + c60.cif2) - 1.0);
  if (unityErr > 1e-12) {
    throw new Error(`Partition of unity violated: ${unityErr}`);
  }

  const diffs = [
    Math.abs(res.shr1 - vec.expected.shr1),
    Math.abs(res.hr2 - vec.expected.hr2),
    Math.abs(c36.os - vec.expected.os_36m),
    Math.abs(c36.cif1 - vec.expected.cif1_36m),
    Math.abs(c60.os - vec.expected.os_60m),
    Math.abs(c60.cif1 - vec.expected.cif1_60m),
    Math.abs(c60.cif2 - vec.expected.cif2_60m),
    Math.abs(c120.os - vec.expected.os_120m)
  ];
  maxDiff = Math.max(maxDiff, ...diffs);
}

console.log(`[Node.js <-> Python Parity Test] Evaluated ${vectors.length} profiles. Max absolute difference = ${maxDiff.toExponential(3)}`);
if (maxDiff > 1e-7) {
  throw new Error(`Parity check failed: maxDiff ${maxDiff} > 1e-7`);
}
console.log('PASS: JavaScript Royston-Parmar + Fine-Gray engine matches Python reference to < 1e-7!');
