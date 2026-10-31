// Generates tests/driftlet_randles.json for tests/test_driftlet_randles.py: spectra simulated by
// driftlet (https://github.com/markblundeberg/driftlet, 0BSD) for its Randles page
// (demos/randles.html: Pt in 1 mM Fe3+/Fe2+ with KCl, at rest, double layer below the grid), and the
// Randles circuit's elements from the page's physics. The physics is imported from the page itself,
// as driftlet's own test/demos.test.js does, so the data are the page's checked theory.
// Needs Node >= 22 and driftlet >= 1.4.0 (before it an electrode's kinetics depended on the mesh).
// Usage: node tests/driftlet_randles.mjs <driftlet checkout> tests/driftlet_randles.json
import { readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { execFileSync } from 'node:child_process';

const [root, out] = process.argv.slice(2);
// The page-free head of the page's script, as driftlet's test/demos.test.js imports it.
const html = readFileSync(join(root, 'demos/randles.html'), 'utf8');
const MARKER = '// (All above runs without a page';
const script = /<script type="module">([\s\S]*?)<\/script>/.exec(html)?.[1];
if (!script || !script.includes(MARKER)) throw new Error(`demos/randles.html: no module script with the marker "${MARKER}"`);
const head = script
  .slice(0, script.indexOf(MARKER))
  .replace(/^import .*'\.\/(plot|common)\.js';$/gm, '')
  .replaceAll("'../src/", `'${pathToFileURL(join(root, 'src/')).href}`);
if (/document\./.test(head)) throw new Error('demos/randles.html: page code above the marker');
// T (K) and c (mol/m3 of each iron ion) are the page's own, so the per-ion Warburgs below use its values.
const { atRest, randles, T, c } = await import(`data:text/javascript;base64,${Buffer.from(`${head}\nexport { atRest, randles, T, c };`).toString('base64')}`);
const { IONS } = await import(pathToFileURL(join(root, 'src/kit.js')).href);
const { GAS_CONSTANT, FARADAY } = await import(pathToFileURL(join(root, 'src/index.js')).href);

const RT = GAS_CONSTANT * T;
const fs = Array.from({ length: 81 }, (_, k) => 10 ** (-2 + k / 10)); // 10 mHz to 1 MHz, 10 per decade
// The page test's corners and middle, plus the page's default (k0 cm/s, salt mol/m3, L m)
const sets = [
  { k0: 1e-2, salt: 100, L: 10 ** 1.75 * 1e-6 },
  { k0: 1, salt: 10, L: 3.2e-6 },
  { k0: 1, salt: 10, L: 316e-6 },
  { k0: 1e-4, salt: 1000, L: 316e-6 },
  { k0: 1e-4, salt: 10, L: 3.2e-6 },
  { k0: 1e-2, salt: 100, L: 17.8e-6 },
];
const spectra = sets.map((p) => {
  const dev = atRest(p), th = randles(p), r = dev.impedance(fs);
  const z = fs.map(th.Z);
  // Each ion's finite-length Warburg as our Ws: R_W = RT L / (F^2 c D), tau_W = L^2 / D
  const ws = (ion) => ({ R_W: (RT / (FARADAY * FARADAY * c)) * p.L / IONS[ion].D, tau_W: p.L * p.L / IONS[ion].D });
  return {
    params: p,
    physics: { Rs: th.Rs, Rct: th.Rct, Cdl: th.Cdl, RW: th.RW, Ws_Fe3: ws('Fe3+'), Ws_Fe2: ws('Fe2+') },
    Z_re: Array.from(r.Z.re), Z_im: Array.from(r.Z.im),
    circuit_re: z.map((w) => w.re), circuit_im: z.map((w) => w.im),
  };
});
const pkg = JSON.parse(readFileSync(join(root, 'package.json'), 'utf8'));
let source = pkg.version;
try {
  source += ` (${execFileSync('git', ['-C', root, 'rev-parse', '--short', 'HEAD'], { stdio: ['ignore', 'pipe', 'ignore'] }).toString().trim()})`;
} catch {
  // not a git checkout (an npm install, a tarball): the version alone labels the data
}
writeFileSync(out, JSON.stringify({ driftlet: source, units: 'Ohm m^2, F/m^2, s', freq: fs, spectra }));
