/* ==========================================================================
 * Bobo Plan2BOQ — Calculation Engine (pure logic, no DOM)
 * Pipeline: Drawing → Walls/Openings → Rooms → RC Frame (prelim) → QTO → BOQ
 * Units: geometry in m; loads kg/m², forces kgf; member sizes cm; fc', fy in ksc
 * Codes: ACI 318-19 (min depth Table 9.3.1.1 / 7.3.1.1 / 8.3.1.2, φ factors,
 *        punching 22.6.5.2) — preliminary sizing only, verify by full analysis.
 * ========================================================================== */
const Engine = (() => {
  'use strict';

  /* ---------------------------------------------------------------- 1. RASTER */
  function toGray(img) {
    const n = img.w * img.h, g = new Uint8Array(n), d = img.data;
    for (let i = 0; i < n; i++) g[i] = (d[4 * i] * 0.299 + d[4 * i + 1] * 0.587 + d[4 * i + 2] * 0.114) | 0;
    return g;
  }
  function otsu(g) {
    const h = new Float64Array(256); for (const v of g) h[v]++;
    const n = g.length; let sum = 0; for (let i = 0; i < 256; i++) sum += i * h[i];
    let sB = 0, wB = 0, best = 0, th = 128;
    for (let t = 0; t < 256; t++) {
      wB += h[t]; if (!wB) continue; const wF = n - wB; if (!wF) break;
      sB += t * h[t]; const mB = sB / wB, mF = (sum - sB) / wF, v = wB * wF * (mB - mF) ** 2;
      if (v > best) { best = v; th = t; }
    }
    return th;
  }
  /* separable box filter on binary mask: erode = all-ones in window, dilate = any */
  function box1D(src, w, h, k, horiz, erode) {
    const out = new Uint8Array(w * h); const a = (k - 1) >> 1, b = k - 1 - a;
    const lo = erode ? a : b, hi = erode ? b : a; // dilate uses mirrored window
    const L = horiz ? w : h, M = horiz ? h : w; const ps = new Int32Array(L + 1);
    for (let m = 0; m < M; m++) {
      for (let l = 0; l < L; l++) ps[l + 1] = ps[l] + src[horiz ? m * w + l : l * w + m];
      for (let l = 0; l < L; l++) {
        const s = l - lo, e = l + hi; let v;
        if (erode) v = (s >= 0 && e < L) ? (ps[e + 1] - ps[s] === k ? 1 : 0) : 0;
        else { const s2 = Math.max(0, s), e2 = Math.min(L - 1, e); v = ps[e2 + 1] - ps[s2] > 0 ? 1 : 0; }
        out[horiz ? m * w + l : l * w + m] = v;
      }
    }
    return out;
  }
  function opening(mask, w, h, kx, ky) {
    let m = mask;
    if (kx > 1) m = box1D(m, w, h, kx, true, true);
    if (ky > 1) m = box1D(m, w, h, ky, false, true);
    if (kx > 1) m = box1D(m, w, h, kx, true, false);
    if (ky > 1) m = box1D(m, w, h, ky, false, false);
    return m;
  }
  const count = m => { let s = 0; for (let i = 0; i < m.length; i++) s += m[i]; return s; };

  /* Auto noise-kernel: first plateau of dark-pixel survival curve f(s) */
  function autoKernel(mask, w, h) {
    const f = [count(mask)], curve = [];
    for (let s = 2; s <= 14; s++) f.push(count(opening(mask, w, h, s, s)));
    for (let i = 0; i < f.length - 1; i++) curve.push({ s: i + 1, f: f[i], drop: (f[i] - f[i + 1]) / Math.max(1, f[i]) });
    let k = 3;
    for (let i = 1; i < curve.length; i++) if (curve[i].drop < 0.03) { k = curve[i].s; break; }
    return { k: Math.max(2, k), curve };
  }
  function runStats(mask, w, h, maxRun) {
    const hist = new Float64Array(maxRun + 1);
    const scan = (horiz) => {
      const L = horiz ? w : h, M = horiz ? h : w;
      for (let m = 0; m < M; m++) { let r = 0;
        for (let l = 0; l <= L; l++) {
          const v = l < L ? mask[horiz ? m * w + l : l * w + m] : 0;
          if (v) r++; else { if (r && r <= maxRun) hist[r] += r; r = 0; }
        } }
    };
    scan(true); scan(false);
    const tot = hist.reduce((a, b) => a + b, 0); let acc = 0, p50 = 0, p90 = 0;
    for (let r = 1; r <= maxRun; r++) { acc += hist[r]; if (!p50 && acc >= 0.5 * tot) p50 = r; if (!p90 && acc >= 0.9 * tot) { p90 = r; break; } }
    return { p50, p90 };
  }
  function components(mask, w, h) {
    const lab = new Int32Array(w * h), comps = []; const st = new Int32Array(w * h);
    let id = 0;
    for (let i = 0; i < w * h; i++) {
      if (!mask[i] || lab[i]) continue; id++; let sp = 0; st[sp++] = i; lab[i] = id; const px = [];
      while (sp) { const p = st[--sp]; px.push(p); const x = p % w, y = (p / w) | 0;
        if (x > 0 && mask[p - 1] && !lab[p - 1]) { lab[p - 1] = id; st[sp++] = p - 1; }
        if (x < w - 1 && mask[p + 1] && !lab[p + 1]) { lab[p + 1] = id; st[sp++] = p + 1; }
        if (y > 0 && mask[p - w] && !lab[p - w]) { lab[p - w] = id; st[sp++] = p - w; }
        if (y < h - 1 && mask[p + w] && !lab[p + w]) { lab[p + w] = id; st[sp++] = p + w; } }
      comps.push(px);
    }
    return comps;
  }
  /* component → centerline segments by weighted LSQ (handles slight skew of photos).
     Splits where stroke thickness jumps ×2 (thick wall continuing into thin railing line). */
  function compToSegs(px, w, horiz) {
    const acc = new Map();
    for (const p of px) { const x = p % w, y = (p / w) | 0; const u = horiz ? x : y, v = horiz ? y : x;
      const e = acc.get(u) || [0, 0]; e[0]++; e[1] += v; acc.set(u, e); }
    const cols = [...acc.entries()].sort((a, b) => a[0] - b[0]); const groups = []; let gcur = [];
    for (const c of cols) { if (gcur.length) { const ref = gcur.slice(-12).map(q => q[1][0]).sort((a, b) => a - b)[Math.min(5, gcur.slice(-12).length >> 1)];
        const prevU = gcur[gcur.length - 1][0]; if (c[0] - prevU > 2 || Math.max(c[1][0], ref) > 2.2 * Math.min(c[1][0], ref)) { groups.push(gcur); gcur = []; } }
      gcur.push(c); }
    if (gcur.length) groups.push(gcur);
    return groups.filter(g => g.length >= 3).map(g => {
      let sw = 0, su = 0, sv = 0, suu = 0, suv = 0; const cnts = [];
      for (const [u, [c, s]] of g) { const v = s / c; sw += c; su += c * u; sv += c * v; suu += c * u * u; suv += c * u * v; cnts.push(c); }
      const den = sw * suu - su * su; const b = den ? (sw * suv - su * sv) / den : 0, a = (sv - b * su) / sw;
      cnts.sort((p, q) => p - q); const t = cnts[cnts.length >> 1]; const umin = g[0][0], umax = g[g.length - 1][0] + 1; const fn = u => a + b * u;
      return horiz ? { x1: umin, y1: fn(umin), x2: umax, y2: fn(umax), t } : { x1: fn(umin), y1: umin, x2: fn(umax), y2: umax, t };
    });
  }

  /** analyzeRaster: image → wall centerlines (px). opts.k (noise kernel) optional */
  function analyzeRaster(img, opts = {}) {
    const { w, h } = img; const g = toGray(img); const th = opts.threshold || otsu(g);
    // adaptive threshold (local mean, integral image) → robust to photo shading; global guard for flat CAD scans
    const R = Math.max(15, Math.round(Math.max(w, h) / 24)), C = opts.C || 18; const I = new Float64Array((w + 1) * (h + 1));
    for (let y = 0; y < h; y++) { let rs = 0; for (let x = 0; x < w; x++) { rs += g[y * w + x]; I[(y + 1) * (w + 1) + x + 1] = I[y * (w + 1) + x + 1] + rs; } }
    const dark = new Uint8Array(w * h);
    for (let y = 0; y < h; y++) { const y0 = Math.max(0, y - R), y1 = Math.min(h, y + R + 1);
      for (let x = 0; x < w; x++) { const x0 = Math.max(0, x - R), x1 = Math.min(w, x + R + 1);
        const m = (I[y1 * (w + 1) + x1] - I[y0 * (w + 1) + x1] - I[y1 * (w + 1) + x0] + I[y0 * (w + 1) + x0]) / ((x1 - x0) * (y1 - y0));
        const v = g[y * w + x]; dark[y * w + x] = (v < m - C && v < Math.max(th, 200)) || v < th * 0.5 ? 1 : 0; } }
    const ak = autoKernel(dark, w, h); const k = opts.k || ak.k;
    const clean = opening(dark, w, h, k, k);
    const rs = runStats(clean, w, h, 120);
    const L = opts.minLen || Math.max(20, Math.round(1.4 * rs.p50));
    const Hm = opening(clean, w, h, L, 1), Vm = opening(clean, w, h, 1, L);
    const segs = [];
    for (const [m, hz] of [[Hm, true], [Vm, false]]) for (const c of components(m, w, h)) {
      if (c.length < L * 2) continue; for (const s of compToSegs(c, w, hz)) {
      const len = Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
      if (len < L || s.t > 0.75 * len || s.t > 3 * rs.p90 + 4) continue;
      if (s.t < Math.max(opts.minT || 4, 0.25 * rs.p50)) continue; // railing / hatch / furniture line
      segs.push(s); }
    }
    return { w, h, threshold: th, k, kCurve: ak.curve, tP50: rs.p50, tP90: rs.p90, minLen: L, segs, dark, clean };
  }
  /** ratio of dark pixels of ORIGINAL mask inside an opening band → window glazing lines */
  function bandDarkRatio(dark, w, h, op, sPx) {
    const ux = (op.x2 - op.x1), uy = (op.y2 - op.y1), len = Math.hypot(ux, uy) || 1; const nx = -uy / len, ny = ux / len;
    const half = op.t / 2 / sPx; let n = 0, d = 0;
    for (let a = 0.15; a <= 0.85; a += 0.02) for (let b = -half; b <= half; b += 0.7) {
      const x = Math.round((op.x1 + ux * a) / sPx + nx * b), y = Math.round((op.y1 + uy * a) / sPx + ny * b);
      if (x < 0 || y < 0 || x >= w || y >= h) continue; n++; d += dark[y * w + x]; }
    return n ? d / n : 0;
  }

  /* ---------------------------------------------------------------- 2. DXF */
  function parseDXF(text) {
    const L = text.split(/\r?\n/); const pairs = [];
    for (let i = 0; i + 1 < L.length; i += 2) pairs.push([parseInt(L[i].trim(), 10), L[i + 1].trim()]);
    let insunits = 0; for (let i = 0; i < pairs.length - 1; i++) if (pairs[i][1] === '$INSUNITS') { insunits = parseInt(pairs[i + 1][1], 10); break; }
    let s = pairs.findIndex(p => p[0] === 2 && p[1] === 'ENTITIES'); const ents = [];
    if (s < 0) return { insunits, entities: ents, layers: [] };
    let cur = null;
    for (let i = s + 1; i < pairs.length; i++) { const [c, v] = pairs[i];
      if (c === 0) { if (cur) ents.push(cur); if (v === 'ENDSEC') { cur = null; break; } cur = { type: v, layer: '0', pts: [], closed: false }; continue; }
      if (!cur) continue;
      if (c === 8) cur.layer = v;
      else if (c === 10) cur.pts.push([+v, 0]); else if (c === 20 && cur.pts.length) cur.pts[cur.pts.length - 1][1] = +v;
      else if (c === 11) cur.p2 = [+v, 0]; else if (c === 21 && cur.p2) cur.p2[1] = +v;
      else if (c === 40) cur.r = +v; else if (c === 70) cur.closed = (+v & 1) === 1;
    }
    const layers = [...new Set(ents.map(e => e.layer))];
    return { insunits, entities: ents, layers };
  }
  const UNIT_M = { 1: 0.0254, 2: 0.3048, 4: 0.001, 5: 0.01, 6: 1 };
  function dxfSegments(dxf, layerFilter) {
    const segs = [];
    for (const e of dxf.entities) { if (layerFilter && !layerFilter(e.layer)) continue;
      if (e.type === 'LINE' && e.pts[0] && e.p2) segs.push({ x1: e.pts[0][0], y1: e.pts[0][1], x2: e.p2[0], y2: e.p2[1], layer: e.layer });
      if (e.type === 'LWPOLYLINE') { const p = e.pts; for (let i = 0; i < p.length - 1; i++) segs.push({ x1: p[i][0], y1: p[i][1], x2: p[i + 1][0], y2: p[i + 1][1], layer: e.layer });
        if (e.closed && p.length > 2) segs.push({ x1: p[p.length - 1][0], y1: p[p.length - 1][1], x2: p[0][0], y2: p[0][1], layer: e.layer }); }
    }
    return segs;
  }
  function dxfUnitScale(dxf, segs, choice) {
    if (choice && choice !== 'auto') return { s: { mm: 0.001, cm: 0.01, m: 1, ft: 0.3048, in: 0.0254 }[choice], note: 'ผู้ใช้กำหนดหน่วย ' + choice };
    let ext = 0; for (const s of segs) ext = Math.max(ext, Math.abs(s.x2 - s.x1), Math.abs(s.y2 - s.y1));
    const u = UNIT_M[dxf.insunits];
    if (u && ext * u >= 3 && ext * u <= 500) return { s: u, note: '$INSUNITS=' + dxf.insunits };
    for (const [k, v] of [['mm', 0.001], ['cm', 0.01], ['m', 1]]) if (ext * v >= 3 && ext * v <= 500) return { s: v, note: `$INSUNITS=${dxf.insunits} ไม่สมเหตุผล (ขนาด ${ext.toFixed(1)}) → สมมติหน่วย ${k}` };
    return { s: 1, note: 'สมมติหน่วย m' };
  }

  /* ---------------------------------------------------------------- 3. GEOMETRY */
  const len = s => Math.hypot(s.x2 - s.x1, s.y2 - s.y1);
  const ang = s => { let a = Math.atan2(s.y2 - s.y1, s.x2 - s.x1); if (a < 0) a += Math.PI; if (a >= Math.PI - 1e-9) a -= Math.PI; return a; };
  const angDiff = (a, b) => { const d = Math.abs(a - b) % Math.PI; return Math.min(d, Math.PI - d); };
  function lineX(s, q) { // intersection of infinite lines
    const d1x = s.x2 - s.x1, d1y = s.y2 - s.y1, d2x = q.x2 - q.x1, d2y = q.y2 - q.y1; const den = d1x * d2y - d1y * d2x;
    if (Math.abs(den) < 1e-12) return null; const t = ((q.x1 - s.x1) * d2y - (q.y1 - s.y1) * d2x) / den;
    return { x: s.x1 + d1x * t, y: s.y1 + d1y * t };
  }
  function projT(s, p) { const dx = s.x2 - s.x1, dy = s.y2 - s.y1, L2 = dx * dx + dy * dy || 1; return ((p.x - s.x1) * dx + (p.y - s.y1) * dy) / L2; }
  function distPtSeg(p, s) { const t = Math.max(0, Math.min(1, projT(s, p))); return Math.hypot(s.x1 + (s.x2 - s.x1) * t - p.x, s.y1 + (s.y2 - s.y1) * t - p.y); }

  /** snap endpoints onto centerlines of crossing walls (fixes overlapping corner blocks) */
  function snap(walls) {
    for (let it = 0; it < 2; it++) for (const s of walls) for (const end of [1, 2]) {
      const P = { x: s['x' + end], y: s['y' + end] }; let best = null, bd = 1e9;
      for (const q of walls) { if (q === s || angDiff(ang(s), ang(q)) < 0.5) continue;
        const X = lineX(s, q); if (!X) continue; const tol = Math.max(s.t, q.t) * 1.05 + 0.05;
        const d = Math.hypot(X.x - P.x, X.y - P.y); if (d > tol) continue;
        const tq = projT(q, X), Lq = len(q); if (tq * Lq < -tol || (tq - 1) * Lq > tol) continue;
        if (d < bd) { bd = d; best = X; } }
      if (best) { s['x' + end] = best.x; s['y' + end] = best.y; }
    }
    return walls.filter(s => len(s) > 0.05);
  }
  /** merge collinear walls; gaps in [minOpen,maxOpen] become openings */
  function mergeCollinear(walls, o) {
    const lines = [], used = new Set();
    for (let i = 0; i < walls.length; i++) { if (used.has(i)) continue; const base = walls[i]; const grp = [base]; used.add(i);
      for (let j = i + 1; j < walls.length; j++) { if (used.has(j)) continue; const q = walls[j];
        if (angDiff(ang(base), ang(q)) > 0.06) continue;
        const off = Math.max(distLine(base, { x: q.x1, y: q.y1 }), distLine(base, { x: q.x2, y: q.y2 }));
        if (off > Math.max(base.t, q.t) * 0.6 + 0.03) continue; grp.push(q); used.add(j); }
      // project onto base direction
      const L0 = len(base), ux = (base.x2 - base.x1) / L0, uy = (base.y2 - base.y1) / L0;
      const iv = grp.map(q => { const a = (q.x1 - base.x1) * ux + (q.y1 - base.y1) * uy, b = (q.x2 - base.x1) * ux + (q.y2 - base.y1) * uy; return { a: Math.min(a, b), b: Math.max(a, b), t: q.t, src: q }; }).sort((p, q) => p.a - q.a);
      // cluster into chains: gap ≤ maxOpen keeps the same line (beam line continues over openings)
      let chain = [iv[0]];
      const flush = () => { lines.push(makeLine(base, ux, uy, chain, o)); };
      for (let k = 1; k < iv.length; k++) { const lastB = Math.max(...chain.map(c => c.b));
        if (iv[k].a - lastB <= o.maxOpen) chain.push(iv[k]); else { flush(); chain = [iv[k]]; } }
      flush();
    }
    return lines;
  }
  function distLine(s, p) { const L = len(s) || 1; return Math.abs((s.x2 - s.x1) * (s.y1 - p.y) - (s.x1 - p.x) * (s.y2 - s.y1)) / L; }
  function makeLine(base, ux, uy, chain, o) {
    const pt = a => ({ x: base.x1 + ux * a, y: base.y1 + uy * a });
    let pieces = [], openings = []; let cur = { a: chain[0].a, b: chain[0].b, tw: chain[0].t * (chain[0].b - chain[0].a), lw: chain[0].b - chain[0].a };
    for (let k = 1; k < chain.length; k++) { const c = chain[k]; const gap = c.a - cur.b;
      if (gap <= o.minGap) { cur.b = Math.max(cur.b, c.b); cur.tw += c.t * (c.b - c.a); cur.lw += c.b - c.a; }
      else { pieces.push(cur); if (gap >= o.minOpen) openings.push({ a: cur.b, b: c.a }); else pieces.push({ a: cur.b, b: c.a, tw: 0, lw: 0, fill: true });
        cur = { a: c.a, b: c.b, tw: c.t * (c.b - c.a), lw: c.b - c.a }; } }
    pieces.push(cur);
    const t = pieces.reduce((s, p) => s + p.tw, 0) / Math.max(1e-9, pieces.reduce((s, p) => s + p.lw, 0));
    const A = pt(chain[0].a), B = pt(Math.max(...chain.map(c => c.b)));
    return { x1: A.x, y1: A.y, x2: B.x, y2: B.y, t, ux, uy,
      pieces: pieces.map(p => { const P = pt(p.a), Q = pt(p.b); return { x1: P.x, y1: P.y, x2: Q.x, y2: Q.y, t }; }),
      openings: openings.map(g => { const P = pt(g.a), Q = pt(g.b); return { x1: P.x, y1: P.y, x2: Q.x, y2: Q.y, w: g.b - g.a, t }; }) };
  }

  /** remove small isolated groups (furniture, kitchen island, legends) */
  function dropIslands(lines, maxDim, log) {
    const n = lines.length, par = [...Array(n).keys()]; const fd = i => par[i] === i ? i : (par[i] = fd(par[i]));
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) { const a = lines[i], b = lines[j]; const tol = Math.max(a.t, b.t) / 2 + 0.1;
      if (distPtSeg({ x: a.x1, y: a.y1 }, b) < tol || distPtSeg({ x: a.x2, y: a.y2 }, b) < tol || distPtSeg({ x: b.x1, y: b.y1 }, a) < tol || distPtSeg({ x: b.x2, y: b.y2 }, a) < tol) par[fd(i)] = fd(j); }
    const box = {}; lines.forEach((l, i) => { const r = fd(i), b = box[r] || (box[r] = [1e9, 1e9, -1e9, -1e9]); b[0] = Math.min(b[0], l.x1, l.x2); b[1] = Math.min(b[1], l.y1, l.y2); b[2] = Math.max(b[2], l.x1, l.x2); b[3] = Math.max(b[3], l.y1, l.y2); });
    const keep = lines.filter((l, i) => { const b = box[fd(i)]; return Math.max(b[2] - b[0], b[3] - b[1]) >= maxDim; });
    if (keep.length < n) log.push(`ตัดเส้นโดดเดี่ยวขนาดเล็ก (เฟอร์นิเจอร์/เคาน์เตอร์) ${n - keep.length} เส้น`);
    return keep;
  }

  /** free wall ends: extend along axis to the next crossing wall (≤ maxOpen) → gap = opening at corner */
  function extendFreeEnds(lines, P, log) {
    let n = 0;
    for (const l of lines) for (const e of [1, 2]) {
      const p = { x: l['x' + e], y: l['y' + e] }; const tolT = q => Math.max(l.t, q.t) / 2 + 0.08;
      if (lines.some(q => q !== l && distPtSeg(p, q) <= tolT(q))) continue;
      const dir = e === 2 ? 1 : -1; let best = null, bd = 1e9;
      for (const q of lines) { if (q === l || angDiff(ang(l), ang(q)) < 0.5) continue; const X = lineX(l, q); if (!X) continue;
        const d = ((X.x - p.x) * l.ux + (X.y - p.y) * l.uy) * dir; if (d <= 0 || d > P.maxOpen) continue;
        const tq = projT(q, X), Lq = len(q), tl = q.t / 2 + 0.6; if (tq * Lq < -tl || (tq - 1) * Lq > tl) continue;
        if (d < bd) { bd = d; best = X; } }
      if (!best) continue; n++;
      const piece = { x1: p.x, y1: p.y, x2: best.x, y2: best.y, w: bd, t: l.t };
      if (bd - l.t / 2 >= P.minOpen) l.openings.push(piece); else l.pieces.push({ x1: p.x, y1: p.y, x2: best.x, y2: best.y, t: l.t });
      l['x' + e] = best.x; l['y' + e] = best.y;
    }
    if (n) log.push(`ต่อปลายผนังอิสระไปชนผนังตั้งฉาก ${n} จุด (ช่องเปิดชิดมุม)`);
  }

  /* ---------------------------------------------------------------- 4. ROOMS (grid) */
  function buildGrid(lines, cellHint) {
    let xmin = 1e9, ymin = 1e9, xmax = -1e9, ymax = -1e9;
    for (const l of lines) { xmin = Math.min(xmin, l.x1, l.x2); xmax = Math.max(xmax, l.x1, l.x2); ymin = Math.min(ymin, l.y1, l.y2); ymax = Math.max(ymax, l.y1, l.y2); }
    const pad = 1.0; xmin -= pad; ymin -= pad; xmax += pad; ymax += pad;
    if (![xmin, xmax, ymin, ymax].every(Number.isFinite)) throw new Error('พิกัดผนังไม่ถูกต้อง (non-finite)');
    const c = Math.max(cellHint || 0.05, Math.max(xmax - xmin, ymax - ymin) / 700);
    const W = Math.ceil((xmax - xmin) / c), H = Math.ceil((ymax - ymin) / c); const g = new Uint8Array(W * H); // 1 wall, 2 ext, 3+ room
    const draw = (s, t) => { const r = Math.max(t / 2, 0.75 * c);
      const x0 = Math.floor((Math.min(s.x1, s.x2) - r - xmin) / c), x1 = Math.ceil((Math.max(s.x1, s.x2) + r - xmin) / c);
      const y0 = Math.floor((Math.min(s.y1, s.y2) - r - ymin) / c), y1 = Math.ceil((Math.max(s.y1, s.y2) + r - ymin) / c);
      for (let y = Math.max(0, y0); y <= Math.min(H - 1, y1); y++) for (let x = Math.max(0, x0); x <= Math.min(W - 1, x1); x++)
        if (distPtSeg({ x: xmin + (x + 0.5) * c, y: ymin + (y + 0.5) * c }, s) <= r) g[y * W + x] = 1; };
    for (const l of lines) { draw(l, l.t); }
    // flood exterior
    const st = []; const push = i => { if (g[i] === 0) { g[i] = 2; st.push(i); } };
    for (let x = 0; x < W; x++) { push(x); push((H - 1) * W + x); } for (let y = 0; y < H; y++) { push(y * W); push(y * W + W - 1); }
    while (st.length) { const p = st.pop(), x = p % W, y = (p / W) | 0; if (x > 0) push(p - 1); if (x < W - 1) push(p + 1); if (y > 0) push(p - W); if (y < H - 1) push(p + W); }
    // label rooms
    const rooms = []; let id = 3;
    for (let i = 0; i < W * H; i++) { if (g[i] !== 0) continue; const cells = [i]; g[i] = id; let k = 0;
      while (k < cells.length) { const p = cells[k++], x = p % W, y = (p / W) | 0;
        for (const q of [x > 0 ? p - 1 : -1, x < W - 1 ? p + 1 : -1, y > 0 ? p - W : -1, y < H - 1 ? p + W : -1]) if (q >= 0 && g[q] === 0) { g[q] = id; cells.push(q); } }
      let per = 0, bx0 = 1e9, by0 = 1e9, bx1 = -1, by1 = -1;
      for (const p of cells) { const x = p % W, y = (p / W) | 0; bx0 = Math.min(bx0, x); bx1 = Math.max(bx1, x); by0 = Math.min(by0, y); by1 = Math.max(by1, y);
        for (const q of [p - 1, p + 1, p - W, p + W]) if (g[q] === 1 || g[q] === 2) per++; }
      rooms.push({ id, area: cells.length * c * c, perim: per * c, cx: xmin + ((bx0 + bx1) / 2 + 0.5) * c, cy: ymin + ((by0 + by1) / 2 + 0.5) * c,
        bw: (bx1 - bx0 + 1) * c, bh: (by1 - by0 + 1) * c, ncell: cells.length }); id++; if (id > 250) break; }
    let foot = 0; for (let i = 0; i < W * H; i++) if (g[i] !== 2) foot++;
    return { g, W, H, c, xmin, ymin, rooms, footprint: foot * c * c };
  }
  const cellAt = (G, x, y) => { const i = Math.floor((x - G.xmin) / G.c), j = Math.floor((y - G.ymin) / G.c); return (i < 0 || j < 0 || i >= G.W || j >= G.H) ? 2 : G.g[j * G.W + i]; };

  /* ---------------------------------------------------------------- 5. RC FRAME */
  function frame(lines, G, P) {
    const tol = l => l.t / 2 + 0.08; const cols = [];
    const addCol = (x, y, why) => { for (const c of cols) if (Math.hypot(c.x - x, c.y - y) < P.colMerge) { c.why.add(why); return c; } const c = { x, y, why: new Set([why]) }; cols.push(c); return c; };
    for (const l of lines) for (const e of [1, 2]) { const p = { x: l['x' + e], y: l['y' + e] }; let touch = 0;
      for (const q of lines) { if (q === l) continue; if (distPtSeg(p, q) <= Math.max(tol(l), tol(q))) touch++; }
      if (touch) addCol(p.x, p.y, touch > 1 ? 'X' : 'L/T'); else if (len(l) >= P.freeEndMin) addCol(p.x, p.y, 'end'); }
    // intermediate columns: span > Lmax
    for (const l of lines) { const on = cols.filter(c => distPtSeg(c, l) <= tol(l) + 0.05).map(c => projT(l, c) * len(l)).sort((a, b) => a - b);
      for (let i = 0; i < on.length - 1; i++) { const sp = on[i + 1] - on[i]; if (sp > P.Lmax) { const n = Math.ceil(sp / P.Lmax);
        for (let k = 1; k < n; k++) { const a = on[i] + sp * k / n; addCol(l.x1 + l.ux * a, l.y1 + l.uy * a, 'span'); } } } }
    // interior grid for large rooms when there are suspended floors / RC roof (secondary frame ≤ Lmax)
    const gridLines = [];
    if (P.floors > 1 || P.roofRC) for (const r of G.rooms) { if (Math.min(r.bw, r.bh) <= P.Lmax * 1.15) continue;
      const x0 = r.cx - r.bw / 2, y0 = r.cy - r.bh / 2, nx = Math.ceil(r.bw / P.Lmax), ny = Math.ceil(r.bh / P.Lmax);
      const inRoom = (x, y) => cellAt(G, x, y) === r.id;
      const xs = [...Array(nx - 1)].map((_, i) => x0 + r.bw * (i + 1) / nx), ys = [...Array(ny - 1)].map((_, i) => y0 + r.bh * (i + 1) / ny);
      for (const x of xs) for (const y of ys) if (inRoom(x, y)) addCol(x, y, 'grid');
      const runs = (fixed, a0, a1, vert) => { const st = G.c * 2; let s0 = null; for (let a = a0; a <= a1 + 1e-9; a += st) { const ok = vert ? inRoom(fixed, a) : inRoom(a, fixed);
          if (ok && s0 === null) s0 = a; if ((!ok || a + st > a1) && s0 !== null) { const e = ok ? a : a - st; if (e - s0 > 0.5) gridLines.push(vert ? { x1: fixed, y1: s0 - r.bh * 0, x2: fixed, y2: e } : { x1: s0, y1: fixed, x2: e, y2: fixed }); s0 = null; } } };
      for (const x of xs) runs(x, y0, y0 + r.bh, true); for (const y of ys) runs(y, x0, x0 + r.bw, false);
      r.gridSpan = Math.min(r.bw / nx, r.bh / ny); }
    // beams = spans between supports along each line
    const beams = [];
    for (const l of lines) { const L = len(l); const on = cols.filter(c => distPtSeg(c, l) <= tol(l) + 0.05).map(c => Math.max(0, Math.min(L, projT(l, c) * L))).sort((a, b) => a - b);
      const pts = [...new Set([0, ...on, L].map(v => +v.toFixed(3)))].sort((a, b) => a - b);
      for (let i = 0; i < pts.length - 1; i++) { const a = pts[i], b = pts[i + 1]; if (b - a < 0.2) continue;
        const sup = on.some(v => Math.abs(v - a) < 0.05) && on.some(v => Math.abs(v - b) < 0.05);
        beams.push({ x1: l.x1 + l.ux * a, y1: l.y1 + l.uy * a, x2: l.x1 + l.ux * b, y2: l.y1 + l.uy * b, L: b - a, cant: !sup }); } }
    for (const gl of gridLines) { const L = len(gl), ux = (gl.x2 - gl.x1) / L, uy = (gl.y2 - gl.y1) / L;
      const on = cols.filter(c => distPtSeg(c, gl) < 0.1).map(c => projT(gl, c) * L).sort((a, b) => a - b); const pts = [0, ...on, L];
      for (let i = 0; i < pts.length - 1; i++) if (pts[i + 1] - pts[i] > 0.2) beams.push({ x1: gl.x1 + ux * pts[i], y1: gl.y1 + uy * pts[i], x2: gl.x1 + ux * pts[i + 1], y2: gl.y1 + uy * pts[i + 1], L: pts[i + 1] - pts[i], cant: false, grid: true }); }
    // tributary area: nearest column (Voronoi) on footprint cells
    for (const c of cols) c.At = 0;
    if (cols.length) for (let j = 0; j < G.H; j++) for (let i = 0; i < G.W; i++) { const v = G.g[j * G.W + i]; if (v === 2) continue;
      const x = G.xmin + (i + 0.5) * G.c, y = G.ymin + (j + 0.5) * G.c; let bi = 0, bd = 1e18;
      for (let k = 0; k < cols.length; k++) { const d = (cols[k].x - x) ** 2 + (cols[k].y - y) ** 2; if (d < bd) { bd = d; bi = k; } }
      cols[bi].At += G.c * G.c; }
    return { cols, beams };
  }
  const ceilTo = (v, s) => Math.ceil(v / s - 1e-9) * s;
  const f = (v, d = 2) => Number(v).toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });

  /** Preliminary design with exhaustive substitution (ACI 318-19, kgf–cm–ksc) */
  function design(fr, G, lines, P) {
    const out = { sheets: [] }; const fc = P.fc, fy = P.fy, fyMPa = fy / 10.197;
    const rooms = G.rooms; const pan = rooms.map(r => r.gridSpan ? { s: r.gridSpan, l: r.gridSpan * 1.0 } : { s: Math.min(r.bw, r.bh), l: Math.max(r.bw, r.bh) });
    const cp0 = pan.length ? pan.reduce((a, b) => b.s > a.s ? b : a) : { s: 4, l: 4 }; const Ls = cp0.s, Ll = cp0.l;
    // --- SLAB
    const beta = Ll / Ls; let hs, sl = [];
    if (P.floors === 1 && !P.roofRC) { hs = P.slabMin; sl.push(`อาคาร 1 ชั้น + หลังคาโครงเบา → ไม่มีพื้นแขวน (มีเฉพาะพื้นบนดิน หนา 10 cm)`, `ใช้ h = ${hs} cm สำหรับคำนวณน้ำหนักเท่านั้น`); }
    else if (beta <= 2) { const hmin = Ls * 100 * (0.8 + fyMPa / 1400) / (36 + 9 * beta); hs = Math.max(P.slabMin, ceilTo(Math.max(hmin, 9), 1));
      sl.push(`Two-way slab (β = Ll/Ls = ${f(Ll)}/${f(Ls)} = ${f(beta)} ≤ 2) — ACI 318-19 Eq. 8.3.1.2(d) [αfm > 2.0]`,
        `h_min = ln(0.8 + fy/1400)/(36 + 9β) = ${f(Ls * 100, 1)}×(0.8 + ${f(fyMPa, 0)}/1400)/(36 + 9×${f(beta)}) = ${f(hmin)} cm ≥ 9 cm`);
    } else { const k = 0.4 + fyMPa / 700; const hmin = Ls * 100 / 24 * k; hs = Math.max(P.slabMin, ceilTo(hmin, 1));
      sl.push(`One-way slab (β = ${f(beta)} > 2) — ACI 318-19 Table 7.3.1.1 (one end continuous) ℓ/24 × (0.4 + fy/700)`, `h_min = ${f(Ls * 100, 1)}/24 × (0.4 + ${f(fyMPa, 0)}/700) = ${f(hmin)} cm`); }
    if (!(P.floors === 1 && !P.roofRC)) sl.push(`ใช้ h = ${hs} cm (ขั้นต่ำตามแนวปฏิบัติ ${P.slabMin} cm) → ${hs >= P.slabMin ? 'PASS ✔' : 'FAIL ✘'}`);
    if (rooms.some(r => r.gridSpan)) sl.push(`ห้องที่ช่วงสั้น > ${P.Lmax} m: เพิ่มเสา/คานภายในตามกริด ≤ Lmax → ช่วงแผ่นพื้นออกแบบ = ${f(Ls)} m`);
    const DLs = hs / 100 * 2400 + P.SDL + P.partition, wu = 1.2 * DLs + 1.6 * P.LL;
    sl.push(`DL = ${f(hs / 100)}×2400 + SDL ${P.SDL} + partition ${P.partition} = ${f(DLs, 0)} kg/m²; LL = ${P.LL} kg/m²`,
      `wu = 1.2DL + 1.6LL = 1.2×${f(DLs, 0)} + 1.6×${P.LL} = ${f(wu, 0)} kg/m²  (ACI 318-19 Eq. 5.3.1b)`);
    out.slab = { h: hs, DL: DLs, wu, Ls, Ll }; out.sheets.push({ title: 'S1 — พื้น (Slab) ห้องที่ช่วงสั้นมากที่สุด', lines: sl });
    // --- BEAM (critical span)
    const sup = fr.beams.filter(b => !b.cant); const bm = sup.length ? sup.reduce((a, b) => b.L > a.L ? b : a) : { L: 4 }; const L = bm.L;
    const hmin = L * 100 / 18.5 * (0.4 + fyMPa / 700); let h = Math.max(P.beamMinH, ceilTo(hmin, 5)); let b = h <= 50 ? 20 : h <= 70 ? 25 : 30; let bl = [];
    bl.push(`ช่วงคานวิกฤต L = ${f(L)} m (คานต่อเนื่องปลายเดียว) — ACI 318-19 Table 9.3.1.1: h_min = ℓ/18.5 × (0.4 + fy/700)`,
      `h_min = ${f(L * 100, 1)}/18.5 × (0.4 + ${f(fyMPa, 0)}/700) = ${f(hmin)} cm`);
    const trib = Math.min(Ls, P.Lmax); let wallH;
    let Mu, phiMn, ok, iter = 0;
    do { wallH = P.Hs - h / 100; const d = h - 6, sw = b / 100 * h / 100 * 2400; const wl = 1.2 * (DLs * trib + sw + P.wallW * wallH) + 1.6 * P.LL * trib;
      Mu = wl * L * L / 10; const rho = P.rhoBeam, As = rho * b * d, a = As * fy / (0.85 * fc * b); phiMn = 0.9 * As * fy * (d - a / 2) / 100;
      ok = phiMn >= Mu;
      if (ok || iter > 12) { bl.push(`ทดลอง b×h = ${b}×${h} cm, d = h − 6 = ${d} cm`,
        `w_u = 1.2(DL·trib + sw + wall) + 1.6(LL·trib) = 1.2(${f(DLs, 0)}×${f(trib)} + ${f(sw, 0)} + ${P.wallW}×${f(wallH)}) + 1.6(${P.LL}×${f(trib)}) = ${f(wl, 0)} kg/m`,
        `M_u = w_u L²/10 = ${f(wl, 0)}×${f(L)}²/10 = ${f(Mu, 0)} kg·m`,
        `As (ρ=${P.rhoBeam}) = ${P.rhoBeam}×${b}×${d} = ${f(As)} cm²; a = As·fy/(0.85fc'b) = ${f(As)}×${fy}/(0.85×${fc}×${b}) = ${f(a)} cm`,
        `φM_n = 0.9×As×fy×(d − a/2) = 0.9×${f(As)}×${fy}×(${d} − ${f(a / 2)})/100 = ${f(phiMn, 0)} kg·m ${ok ? '≥' : '<'} M_u = ${f(Mu, 0)} → ${ok ? 'PASS ✔' : 'FAIL ✘'}`); break; }
      h += 5; b = h <= 50 ? 20 : h <= 70 ? 25 : 30; iter++; } while (true);
    bl.push(`ตรวจความลึกขั้นต่ำ: h = ${h} cm ≥ h_min = ${f(hmin)} cm → ${h >= hmin ? 'PASS ✔' : 'FAIL ✘'}`);
    out.beam = { b, h, L }; out.sheets.push({ title: 'B1 — คาน (Beam) ช่วงยาวที่สุด', lines: bl });
    // --- COLUMN (max tributary)
    const nSus = P.floors - 1; const roofDL = P.roofRC ? DLs : P.roofDL, roofLL = P.roofRC ? P.roofLLrc : P.roofLL;
    const wuRoof = 1.2 * roofDL + 1.6 * roofLL; const crit = fr.cols.length ? fr.cols.reduce((a, c) => c.At > a.At ? c : a) : { At: 10 };
    const At = crit.At; const Lbw = Math.min(2 * Math.sqrt(At), 2 * P.Lmax); const beamW = 1.2 * (b / 100 * h / 100 * 2400) * Lbw; // allowance: beam length per level
    const Pu = At * (nSus * wu + wuRoof) + (nSus + 1) * beamW; const rho = P.rhoCol;
    const k = 0.65 * 0.80 * (0.85 * fc * (1 - rho) + fy * rho); const Ag = Pu / k; let bc = Math.max(P.colMin, ceilTo(Math.sqrt(Ag), 5), b);
    const Puc = Pu + 1.2 * (bc / 100) ** 2 * 2400 * P.Hs * P.floors; let Ag2 = Puc / k; bc = Math.max(bc, ceilTo(Math.sqrt(Ag2), 5)); const phiPn = k * bc * bc;
    out.col = { b: bc, Pu: Puc, At };
    out.sheets.push({ title: 'C1 — เสา (Column) พื้นที่รับน้ำหนักมากที่สุด', lines: [
      `A_t (Voronoi ตามตำแหน่งเสา) = ${f(At)} m²; จำนวนชั้นพื้นแขวน = ${nSus}, หลังคา ${P.roofRC ? 'พื้น คสล.' : 'โครงหลังคาเบา'}`,
      `w_u,roof = 1.2×${f(roofDL, 0)} + 1.6×${roofLL} = ${f(wuRoof, 0)} kg/m²`,
      `น้ำหนักคาน (ยาว min(2√A_t, 2Lmax) = ${f(Lbw)} m ต่อชั้น) = 1.2×${b / 100}×${h / 100}×2400×${f(Lbw)} = ${f(beamW, 0)} kg/ชั้น`,
      `P_u = A_t(n·w_u + w_u,roof) + (n+1)·beam = ${f(At)}(${nSus}×${f(wu, 0)} + ${f(wuRoof, 0)}) + ${nSus + 1}×${f(beamW, 0)} = ${f(Pu, 0)} kgf`,
      `+ น้ำหนักเสา 1.2×(${bc}/100)²×2400×${P.Hs}×${P.floors} → P_u = ${f(Puc, 0)} kgf`,
      `Tied column ACI 318-19 Eq. 22.4.2.2 + Table 21.2.2: φP_n,max = 0.65×0.80[0.85fc'(Ag − Ast) + fy·Ast]`,
      `ρ = ${rho}: φP_n/Ag = 0.52×[0.85×${fc}×(1 − ${rho}) + ${fy}×${rho}] = ${f(k)} ksc`,
      `Ag,req = ${f(Puc, 0)}/${f(k)} = ${f(Ag2, 0)} cm² → ใช้ ${bc}×${bc} cm (Ag = ${bc * bc} cm², ≥ ${P.colMin}×${P.colMin} และ ≥ b คาน)`,
      `φP_n = ${f(k)}×${bc * bc} = ${f(phiPn, 0)} kgf ${phiPn >= Puc ? '≥' : '<'} P_u = ${f(Puc, 0)} → ${phiPn >= Puc ? 'PASS ✔' : 'FAIL ✘'}`] });
    // --- FOOTING per column (service load)
    const wsF = DLs + P.LL, wsR = roofDL + roofLL; const foots = []; let fl = null;
    for (const c of fr.cols) { const Ps = c.At * (nSus * wsF + wsR) + (nSus + 1) * beamW / 1.2 + (bc / 100) ** 2 * 2400 * (P.Hs * P.floors + P.Df);
      const Areq = 1.10 * Ps / (P.qa * 1000); let B = Math.max(P.fMin, ceilTo(Math.sqrt(Areq), 0.10)); let t = 0.30, okp = false, Vu, phiVc, bo, d;
      const Pu_f = 1.45 * Ps; // average load factor (1.2D+1.6L) ≈ 1.45 for D/L mix
      while (!okp && t < 1.5) { d = t * 100 - 7.5; const q = Pu_f / (B * B); const cp = (bc + d) / 100; bo = 4 * (bc + d);
        Vu = q * (B * B - cp * cp); phiVc = 0.75 * 1.06 * Math.sqrt(fc) * bo * d; okp = phiVc >= Vu; if (!okp) t += 0.05; }
      c.Ps = Ps; c.B = B; c.tf = t; foots.push({ B, t });
      if (!fl || Ps > fl.Ps) fl = { Ps, Areq, B, t, Pu_f, d, bo, Vu, phiVc, At: c.At }; }
    if (fl) out.sheets.push({ title: 'F1 — ฐานราก (Isolated Footing) เสาวิกฤต', lines: [
      `P_s = A_t(n(DL+LL) + roof) + beam + เสา = ${f(fl.At)}(${nSus}×${f(wsF, 0)} + ${f(wsR, 0)}) + ... = ${f(fl.Ps, 0)} kgf`,
      `A_req = 1.10 P_s/q_a = 1.10×${f(fl.Ps, 0)}/(${P.qa}×1000) = ${f(fl.Areq, 3)} m² → B×B = ${f(fl.B)}×${f(fl.B)} m (≥ ${P.fMin} m)`,
      `P_u ≈ 1.45 P_s = ${f(fl.Pu_f, 0)} kgf; t = ${f(fl.t)} m, d = t − 7.5 cm = ${f(fl.d, 1)} cm, b_o = 4(c + d) = ${f(fl.bo, 1)} cm`,
      `Punching (ACI 318-19 Table 22.6.5.2(a)): φV_c = 0.75×1.06√fc'·b_o·d = 0.75×1.06×√${fc}×${f(fl.bo, 1)}×${f(fl.d, 1)} = ${f(fl.phiVc, 0)} kgf`,
      `V_u = q_u(B² − (c+d)²) = ${f(fl.Vu, 0)} kgf → ${fl.phiVc >= fl.Vu ? 'PASS ✔' : 'FAIL ✘'}`,
      `q_a = ${P.qa} t/m² เป็นค่าสมมติ — ต้องยืนยันจากผลเจาะสำรวจดิน / Pile load test`] });
    out.foots = foots;
    return out;
  }

  /* ---------------------------------------------------------------- 6. QTO */
  function qto(model, P) {
    const { lines, G, fr, D } = model; const N = P.floors, Hs = P.Hs; const q = [];
    const add = (cat, key, desc, qty, unit, note = '') => q.push({ cat, key, desc, qty, unit, note });
    const nC = fr.cols.length, bc = D.col.b / 100, hb = D.beam.h / 100, bb = D.beam.b / 100, hs = D.slab.h / 100;
    const beamLen = fr.beams.reduce((s, b) => s + b.L, 0), beamLenW = fr.beams.filter(b => !b.grid).reduce((s, b) => s + b.L, 0); const foot = G.footprint;
    // concrete
    const Vcol = nC * bc * bc * (Hs * N + P.Df - 0.4); const Vgb = beamLenW * bb * hb; const Vbm = beamLen * bb * (hb - hs) * (N - 1 + (P.roofRC ? 1 : 0)) + (P.roofRC ? 0 : beamLen * bb * hb * 0.6);
    const Vslab = foot * hs * (N - 1 + (P.roofRC ? 1 : 0)); const Vsog = foot * 0.10; let Vf = 0, Af = 0, Ex = 0;
    for (const c of fr.cols) { Vf += c.B * c.B * c.tf; Af += 4 * c.B * c.tf; Ex += (c.B + 0.6) ** 2 * P.Df; }
    const VexGB = beamLenW * (bb + 0.4) * (hb + 0.1);
    add('2', 'exc', 'ขุดดินฐานราก/คานคอดิน', Ex + VexGB, 'm³', `ฐาน (B+0.6)²×Df, Df=${P.Df} m`);
    add('2', 'fill', 'ถมดินกลับบดอัด', Math.max(0, Ex + VexGB - Vf - Vgb - nC * bc * bc * P.Df), 'm³');
    add('2', 'sand', 'ทรายหยาบรองพื้น หนา 5 cm', fr.cols.reduce((s, c) => s + c.B * c.B, 0) * 0.05 + beamLenW * bb * 0.05, 'm³');
    add('2', 'lean', "คอนกรีตหยาบ (lean) หนา 5 cm", fr.cols.reduce((s, c) => s + c.B * c.B, 0) * 0.05 + beamLenW * bb * 0.05, 'm³');
    add('2', 'c_f', `คอนกรีตฐานราก fc' ${P.fc} ksc`, Vf, 'm³', `${nC} ฐาน`);
    add('2', 'c_c', `คอนกรีตเสา ${D.col.b}×${D.col.b} cm`, Vcol, 'm³', `${nC} ต้น × ${f(Hs * N + P.Df - 0.4)} m`);
    add('2', 'c_gb', `คอนกรีตคานคอดิน ${D.beam.b}×${D.beam.h} cm`, Vgb, 'm³', `ยาวรวม ${f(beamLenW)} m`);
    add('2', 'c_b', `คอนกรีตคานชั้นบน/คานหลังคา`, Vbm, 'm³', P.roofRC ? '' : 'คานรับโครงหลังคา 0.6×');
    add('2', 'c_s', `คอนกรีตพื้นแขวน หนา ${D.slab.h} cm`, Vslab, 'm³', `${f(foot)} m² × ${N - 1 + (P.roofRC ? 1 : 0)} ชั้น`);
    add('2', 'c_sog', 'คอนกรีตพื้นบนดิน หนา 10 cm', Vsog, 'm³');
    const R = P.rebar; const kg = Vf * R.f + Vcol * R.c + (Vgb + Vbm) * R.b + Vslab * R.s + Vsog * R.sog;
    add('2', 'rebar', 'เหล็กเสริม SD40/SR24 (รวม 5% เผื่อ)', kg * 1.05, 'kg', `อัตรา kg/m³: ฐาน ${R.f}, เสา ${R.c}, คาน ${R.b}, พื้น ${R.s}, พื้นบนดิน ${R.sog}`);
    const Fw = Af + nC * 4 * bc * (Hs * N + P.Df - 0.4) + beamLenW * 2 * hb + beamLen * (bb + 2 * (hb - hs)) * (N - 1) + foot * (N - 1 + (P.roofRC ? 1 : 0));
    add('2', 'form', 'ไม้แบบ (ฐาน เสา คาน ท้องพื้น)', Fw, 'm²');
    // architecture
    let wallL = 0, door = [], win = []; for (const l of lines) { for (const p of l.pieces) wallL += len(p); for (const o of l.openings) (o.type === 'window' ? win : door).push(o); }
    const hw = Hs - hb; const colDeduct = nC * bc; const netL = Math.max(0, wallL - colDeduct);
    const winA = win.reduce((s, o) => s + o.w * P.winH, 0); const doorA = door.reduce((s, o) => s + o.w * P.doorH, 0);
    const overOpen = door.reduce((s, o) => s + o.w * Math.max(0, hw - P.doorH), 0) + win.reduce((s, o) => s + o.w * Math.max(0, hw - P.winH), 0);
    const masonry = (netL * hw + overOpen) * N;
    add('3', 'wall', 'ผนังก่ออิฐมวลเบา หนา 7.5 cm', masonry, 'm²', `ยาวสุทธิ ${f(netL)} m × ${f(hw)} m × ${N} ชั้น + เหนือ/ใต้ช่องเปิด`);
    add('3', 'plaster', 'ฉาบปูนผนัง 2 ด้าน', masonry * 2, 'm²');
    add('3', 'paint', 'ทาสีผนัง (รองพื้น + จริง 2 เที่ยว)', masonry * 2, 'm²');
    const roomA = G.rooms.reduce((s, r) => s + r.area, 0), roomP = G.rooms.reduce((s, r) => s + r.perim, 0);
    add('3', 'floor', 'ปูกระเบื้องพื้น', roomA * N, 'm²', `${G.rooms.length} ห้อง`);
    add('3', 'ceil', 'ฝ้ายิปซัมฉาบเรียบ + โครงคร่าว', roomA * N, 'm²');
    add('3', 'skirt', 'บัวเชิงผนัง', Math.max(0, roomP - 2 * door.reduce((s, o) => s + o.w, 0)) * N, 'm');
    add('3', 'door', 'ประตู + วงกบ + อุปกรณ์', door.length * N, 'ชุด', door.map(o => f(o.w)).join(', ') + ' m');
    add('3', 'win', 'หน้าต่างอลูมิเนียม + กระจก', winA * N, 'm²', `${win.length * N} ช่อง สูง ${P.winH} m`);
    add('3', 'roof', P.roofRC ? 'กันซึมดาดฟ้า' : 'หลังคาเมทัลชีท + โครงเหล็ก', foot * (P.roofRC ? 1 : 1.25), 'm²', P.roofRC ? '' : 'รวมชายคา/ความชัน ×1.25');
    add('4', 'elec', 'ระบบไฟฟ้าและแสงสว่าง (เหมาต่อพื้นที่)', foot * N, 'm²');
    add('5', 'sani', 'ระบบประปา-สุขาภิบาล (เหมาต่อพื้นที่)', foot * N, 'm²');
    add('1', 'site', 'ปรับพื้นที่ / ทำผังบริเวณ', foot * 1.5, 'm²', 'พื้นที่ footprint × 1.5');
    add('1', 'temp', 'งานชั่วคราว ป้าย รั้ว ห้องเก็บของ', 1, 'เหมา');
    return { items: q, stats: { wallL, netL, doors: door.length, windows: win.length, roomA, foot, beamLen, nC, masonry } };
  }

  /* ---------------------------------------------------------------- 7. BOQ */
  const CATS = { '1': 'งานเตรียมการ', '2': 'งานโครงสร้าง', '3': 'งานสถาปัตยกรรม', '4': 'งานระบบไฟฟ้า', '5': 'งานระบบสุขาภิบาล' };
  function boq(Q, prices) {
    const rows = []; let grand = 0; const sub = {};
    for (const c of Object.keys(CATS)) { const its = Q.items.filter(i => i.cat === c); let s = 0; let n = 0;
      for (const it of its) { const p = prices[it.key] || [0, 0]; const qty = Math.round(it.qty * 100) / 100; const m = qty * p[0], l = qty * p[1];
        rows.push({ cat: c, no: `${c}.${++n}`, desc: it.desc, qty, unit: it.unit, mu: p[0], mt: m, lu: p[1], lt: l, tot: m + l, note: it.note, key: it.key }); s += m + l; }
      sub[c] = s; grand += s; }
    return { rows, sub, grand, CATS };
  }

  /* ---------------------------------------------------------------- PIPELINE */
  /** walls (m, with t) → full model. opts.classify(op) → 'door'|'window' */
  function buildModel(walls, P, classify) {
    const log = [];
    let W = snap(walls.filter(w => [w.x1, w.y1, w.x2, w.y2, w.t].every(Number.isFinite)).map(w => ({ ...w })));
    let lines = mergeCollinear(W, { minGap: P.minGap, minOpen: P.minOpen, maxOpen: P.maxOpen });
    lines = dropIslands(lines, P.islandMax, log);
    extendFreeEnds(lines, P, log);
    log.push(`Walls: ${walls.length} raw → ${lines.length} wall lines, ${lines.reduce((s, l) => s + l.openings.length, 0)} openings`);
    const G = buildGrid(lines, P.cell);
    G.rooms = G.rooms.filter(r => r.area >= P.minRoom);
    for (const l of lines) { const mx = (l.x1 + l.x2) / 2, my = (l.y1 + l.y2) / 2, nx = -l.uy, ny = l.ux, off = l.t / 2 + 2 * G.c;
      l.ext = cellAt(G, mx + nx * off, my + ny * off) === 2 || cellAt(G, mx - nx * off, my - ny * off) === 2;
      for (const o of l.openings) { o.ext = l.ext; o.type = classify ? classify(o, l) : (l.ext ? 'window' : 'door'); } }
    log.push(`Rooms: ${G.rooms.length}, footprint ${f(G.footprint)} m², grid ${G.W}×${G.H} @ ${f(G.c, 3)} m`);
    const fr = frame(lines, G, P); log.push(`RC frame: ${fr.cols.length} columns, ${fr.beams.length} beam spans`);
    const D = design(fr, G, lines, P);
    const model = { lines, G, fr, D, log }; model.Q = qto(model, P); return model;
  }

  const DEFAULTS = { floors: 2, Hs: 3.0, Df: 1.5, fc: 240, fy: 4000, LL: 200, SDL: 100, partition: 100, wallW: 90, roofDL: 60, roofLL: 30, roofLLrc: 100, roofRC: false,
    qa: 10, Lmax: 5.0, colMerge: 0.6, freeEndMin: 2.0, colMin: 20, beamMinH: 40, slabMin: 10, fMin: 1.0, rhoBeam: 0.008, rhoCol: 0.01,
    minGap: 0.15, minOpen: 0.55, islandMax: 2.0, maxOpen: 3.5, minRoom: 1.0, cell: 0.05, doorH: 2.0, winH: 1.1,
    rebar: { f: 70, c: 160, b: 130, s: 90, sog: 40 } };
  const PRICES = { site: [0, 25], temp: [30000, 20000], exc: [0, 150], fill: [180, 120], sand: [550, 180], lean: [1900, 450], c_f: [2350, 450], c_c: [2350, 650], c_gb: [2350, 550], c_b: [2350, 550], c_s: [2350, 450], c_sog: [2250, 400],
    rebar: [27, 5], form: [220, 160], wall: [230, 120], plaster: [70, 110], paint: [55, 40], floor: [380, 180], ceil: [260, 130], skirt: [90, 40], door: [6500, 1200], win: [2600, 400], roof: [650, 200], elec: [450, 250], sani: [350, 200] };

  return { analyzeRaster, bandDarkRatio, parseDXF, dxfSegments, dxfUnitScale, buildModel, boq, DEFAULTS, PRICES, CATS, fmt: f, len };
})();
if (typeof module !== 'undefined') module.exports = Engine;
