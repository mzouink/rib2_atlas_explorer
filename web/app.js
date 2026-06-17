// RNA Atlas Explorer — client-side filtering/ranking + lazy deep view.
// Sources are chosen as checkboxes (header "Source" menu) and MERGED: FOLDS is the
// union of every checked dataset, each fold tagged with f._dsid so the deep view can
// dispatch struct/ext/reactivity/motifs per row.
let FOLDS = [], MOTIF_SET = [], LETTERS = [];
let MOTIFS_BY_DS = {}, PAIRING_BY_DS = {};  // {dsid: {foldId: ...}} — only motif-bearing datasets populate these
let sortOverride = null;  // {key, dir} from header click
const DATASETS = window.DATASETS || [{ id: "ribo2", label: "curated", base: "", ext: "cif", react: true, motifs: true }];
const DSBYID = {}; DATASETS.forEach((d) => { DSBYID[d.id] = d; });
const LOADED = {};        // dsid -> {folds, motifs, pairing} cache (loaded once)
function dsFor(f) { return DSBYID[f._dsid] || DATASETS[0]; }
function prefix(ds) { return ds && ds.base ? ds.base + "/" : ""; }
function activeSources() { return [...document.querySelectorAll(".src:checked")].map((c) => c.value); }

// Data source: "" = same origin (local serve.py). Otherwise a CloudFront URL
// fronting the S3 data, gated by a shared token (prompted once, kept in localStorage).
const DATA_BASE = (window.DATA_BASE || "").replace(/\/$/, "");
const GATED = !!window.GATED;
function token() { return GATED ? (localStorage.getItem("atlas_token") || "") : ""; }
function durl(path) {
  const base = DATA_BASE ? `${DATA_BASE}/${path}` : path;
  if (!GATED) return base;
  const t = token();
  return t ? `${base}${base.includes("?") ? "&" : "?"}t=${encodeURIComponent(t)}` : base;
}
async function getJSON(path) {
  const r = await fetch(durl(path));
  if (!r.ok) { const e = new Error("http " + r.status); e.status = r.status; throw e; }
  return r.json();
}

const $ = (id) => document.getElementById(id);
const num = (x, d = 1) => (x === null || x === undefined || Number.isNaN(x)) ? "" : (+x).toFixed(d);
// escape for HTML text + attribute values (ids/names can contain " < & — e.g. OpenKnot ids)
const esc = (s) => String(s == null ? "" : s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");

// --- persist filter settings across reloads ---
const FKEY = "atlas_filters";
const FIELD_IDS = ["search", "len_min", "len_max", "plddt_min", "clash_max", "tm_max", "tm_has", "novel_only",
  "ov_max", "shape_ok", "agr_min", "cr_min", "bp_min", "req_tert", "req_rare", "pk", "rank_key", "topn", "per_letter", "alt_palette", "color_by"];
const altPalette = () => !!($("alt_palette") && $("alt_palette").checked);
function snapshot() {
  const s = {};
  FIELD_IDS.forEach((id) => { const el = $(id); if (el) s[id] = el.type === "checkbox" ? el.checked : el.value; });
  s.mf = [...document.querySelectorAll(".mf:checked")].map((c) => c.value);
  s.lf = [...document.querySelectorAll(".lf:checked")].map((c) => c.value);
  s.src = activeSources();
  s.sort = sortOverride;
  s.collapsed = [...document.querySelectorAll("#config fieldset")].map((fs, i) => fs.classList.contains("collapsed") ? i : -1).filter((i) => i >= 0);
  s.allcollapsed = $("config").classList.contains("allcollapsed");
  return s;
}
function applyState(s) {
  if (!s) return;
  FIELD_IDS.forEach((id) => { const el = $(id); if (el && id in s) { if (el.type === "checkbox") el.checked = !!s[id]; else el.value = s[id]; } });
  if (s.mf) document.querySelectorAll(".mf").forEach((c) => { c.checked = s.mf.includes(c.value); });
  if (s.lf) document.querySelectorAll(".lf").forEach((c) => { c.checked = s.lf.includes(c.value); });
  sortOverride = s.sort || null;
  if (s.collapsed) { const fs = document.querySelectorAll("#config fieldset"); s.collapsed.forEach((i) => fs[i] && fs[i].classList.add("collapsed")); }
  if (s.allcollapsed) $("config").classList.add("allcollapsed");
}
function saveState() { try { localStorage.setItem(FKEY, JSON.stringify(snapshot())); } catch (e) {} }
function loadState() { try { return JSON.parse(localStorage.getItem(FKEY) || "null"); } catch (e) { return null; } }

function showGate(msg) {
  const g = $("gate");
  g.classList.remove("hidden");
  $("gate-msg").textContent = msg || "";
  const inp = $("gate-input");
  inp.value = "";
  inp.focus();
  const submit = () => {
    const v = inp.value.trim();
    if (!v) return;
    localStorage.setItem("atlas_token", v);
    g.classList.add("hidden");
    boot();
  };
  $("gate-go").onclick = submit;
  inp.onkeydown = (e) => { if (e.key === "Enter") submit(); };
}

async function boot() {
  if (GATED && !token()) return showGate();
  buildSourcePanel();
  wireSourceUI();
  wireStatic();
  // initial source selection: persisted set, else just the first dataset (ribo2)
  const st = loadState();
  const saved = st && st.src;
  document.querySelectorAll(".src").forEach((c) => {
    c.checked = (saved && saved.length) ? saved.includes(c.value) : (c.value === DATASETS[0].id);
  });
  toggleLetterVisibility(); updateSourceBtn();
  await loadSources();
}

function buildSourcePanel() {
  // The motif-bearing dataset (ribo2) gets the nested per-letter checkboxes (#letter_filter).
  $("source-panel").innerHTML = DATASETS.map((d) => {
    const row = `<label class="srcrow"><input type="checkbox" class="src" value="${d.id}">${d.label}</label>`;
    return d.motifs ? row + `<div id="letter_filter" class="letters srcsub"></div>` : row;
  }).join("");
}

function wireSourceUI() {
  $("source-btn").addEventListener("click", (e) => { e.stopPropagation(); $("source-panel").classList.toggle("hidden"); });
  document.addEventListener("click", (e) => { if (!$("sourcewrap").contains(e.target)) $("source-panel").classList.add("hidden"); });
  document.querySelectorAll(".src").forEach((c) =>
    c.addEventListener("change", () => { toggleLetterVisibility(); updateSourceBtn(); saveState(); loadSources(); }));
}

function toggleLetterVisibility() {
  const on = [...document.querySelectorAll(".src")].some((c) => c.checked && DSBYID[c.value] && DSBYID[c.value].motifs);
  const lf = $("letter_filter"); if (lf) lf.style.display = on ? "" : "none";
}
function updateSourceBtn() {
  const n = activeSources().length;
  $("source-btn").innerHTML = `Source${n ? ` (${n})` : ""} &#9662;`;
}

async function ensureLoaded(dsid) {
  if (LOADED[dsid]) return;
  const ds = DSBYID[dsid];
  const folds = await getJSON(prefix(ds) + "data/folds.json");
  folds.forEach((f) => { f._dsid = dsid; });
  let motifs = {}, pairing = {};
  if (ds.motifs) {
    try { motifs = await getJSON(prefix(ds) + "data/motifs.json"); } catch (e) {}
    try { pairing = await getJSON(prefix(ds) + "data/pairing.json"); } catch (e) {}
  }
  LOADED[dsid] = { folds, motifs, pairing };
}

async function loadSources() {
  const active = activeSources();
  try {
    for (const id of active) await ensureLoaded(id);
  } catch (e) {
    if (GATED && e.status === 403) { localStorage.removeItem("atlas_token"); return showGate("Incorrect passcode — try again."); }
    if (GATED) return showGate("Could not load data (" + (e.status || "network") + ").");
    throw e;
  }
  FOLDS = []; MOTIFS_BY_DS = {}; PAIRING_BY_DS = {}; FBYK = {};
  for (const id of active) {
    const L = LOADED[id]; if (!L) continue;
    FOLDS = FOLDS.concat(L.folds);
    MOTIFS_BY_DS[id] = L.motifs; PAIRING_BY_DS[id] = L.pairing;
  }
  const ms = new Set(), ls = new Set();
  let maxLen = 0;
  for (const f of FOLDS) {
    (f.motifs || []).forEach((m) => ms.add(m));
    if (DSBYID[f._dsid] && DSBYID[f._dsid].motifs && f.letter) ls.add(f.letter);  // letters only from motif-bearing source
    if (f.length > maxLen) maxLen = f.length;
  }
  MOTIF_SET = [...ms].sort();
  LETTERS = [...ls].sort();
  $("len_max").value = maxLen || 9999;
  const maxCR = Math.max(1, ...FOLDS.map((f) => f.contact_ratio || 0));
  $("cr_min").max = Math.ceil(maxCR * 20) / 20;
  buildMotifFilter();
  buildLetterFilter();
  wireDynamic();
  applyState(loadState());
  toggleLetterVisibility();
  syncLabels();
  render();
}

function buildMotifFilter() {
  $("motif_filter").innerHTML = MOTIF_SET.map((m) =>
    `<label><input type="checkbox" class="mf" value="${m}">` +
    `<span class="motif-chip" style="background:${motifColor(m)}">${m.replace(/_/g, " ").toLowerCase()}</span></label>`
  ).join("");
}
function buildLetterFilter() {
  $("letter_filter").innerHTML = LETTERS.map((l) =>
    `<label><input type="checkbox" class="lf" value="${l}" checked>${l}</label>`).join("");
}

function wireDynamic() {
  document.querySelectorAll(".mf,.lf").forEach((c) =>
    c.addEventListener("change", () => { saveState(); render(); }));
}

function wireStatic() {
  document.querySelectorAll("#config input:not(.mf):not(.lf), #config select").forEach((el) =>
    el.addEventListener("input", () => { syncLabels(); saveState(); render(); }));
  $("reset").addEventListener("click", () => {
    document.querySelectorAll(".mf").forEach((c) => c.checked = false);
    document.querySelectorAll(".lf").forEach((c) => c.checked = true);
    ["plddt_min", "tm_max", "ov_max"].forEach((id) => $(id).value = $(id).max);
    $("agr_min").value = -1;
    $("plddt_min").value = 0; $("clash_max").value = 9999; $("len_min").value = 0;
    $("len_max").value = Math.max(...FOLDS.map((f) => f.length || 0));
    ["shape_ok", "req_tert", "req_rare", "tm_has", "novel_only", "per_letter"].forEach((id) => $(id).checked = false);
    $("cr_min").value = 0; $("bp_min").value = 0;
    $("pk").value = "any"; $("rank_key").value = "best_tm1:asc"; $("topn").value = 200;
    if ($("color_by")) $("color_by").value = "a23";
    if ($("search")) $("search").value = "";
    sortOverride = null; localStorage.removeItem(FKEY); syncLabels(); render();
  });
  if ($("search")) $("search").addEventListener("input", () => { saveState(); render(); });
  if ($("alt_palette")) $("alt_palette").addEventListener("change", () => { if (currentDeep) { drawTracks(currentDeep.f, currentDeep.react); load3D(currentDeep.f, currentDeep.react); } });
  if ($("color_by")) $("color_by").addEventListener("change", () => { if (currentDeep) load3D(currentDeep.f, currentDeep.react); });
  // collapsible sections (click a legend) + collapse-all (click the panel heading)
  document.querySelectorAll("#config fieldset legend").forEach((lg) =>
    lg.addEventListener("click", () => { lg.parentElement.classList.toggle("collapsed"); saveState(); }));
  const h2 = document.querySelector("#config h2");
  if (h2) h2.addEventListener("click", () => { $("config").classList.toggle("allcollapsed"); saveState(); });
  $("deep-close").addEventListener("click", closeDeep);
  $("deep").addEventListener("click", (e) => { if (e.target.id === "deep") closeDeep(); });
  document.querySelectorAll('#layoutctl button[data-mode]').forEach((b) =>
    b.addEventListener("click", () => setDeepMode(b.dataset.mode)));
  if ($("help-btn")) $("help-btn").addEventListener("click", () => $("help").classList.remove("hidden"));
  if ($("help-close")) $("help-close").addEventListener("click", () => $("help").classList.add("hidden"));
  $("help").addEventListener("click", (e) => { if (e.target.id === "help") $("help").classList.add("hidden"); });
  updateLayout();
  syncLabels();
}
function syncLabels() {
  $("plddt_min_v").textContent = $("plddt_min").value;
  $("tm_max_v").textContent = (+$("tm_max").value).toFixed(2);
  $("ov_max_v").textContent = (+$("ov_max").value).toFixed(2);
  $("agr_min_v").textContent = (+$("agr_min").value).toFixed(2);
  $("cr_min_v").textContent = (+$("cr_min").value).toFixed(2);
  $("bp_min_v").textContent = (+$("bp_min").value).toFixed(2);
}

function filters() {
  const lf = [...document.querySelectorAll(".lf:checked")].map((c) => c.value);
  const mf = [...document.querySelectorAll(".mf:checked")].map((c) => c.value);
  return {
    q: ($("search") ? $("search").value : "").trim().toLowerCase(),
    lmin: +$("len_min").value, lmax: +$("len_max").value,
    plddt: +$("plddt_min").value, clash: +$("clash_max").value,
    tmax: +$("tm_max").value, tmhas: $("tm_has").checked, novelOnly: $("novel_only").checked,
    ovmax: +$("ov_max").value,
    shape: $("shape_ok").checked, agr: +$("agr_min").value,
    crmin: +$("cr_min").value, bpmin: +$("bp_min").value,
    tert: $("req_tert").checked, rare: $("req_rare").checked, motifs: mf,
    pk: $("pk").value, letters: new Set(lf),
    rank: $("rank_key").value, topn: +$("topn").value, perLetter: $("per_letter").checked,
  };
}

function pass(f, c) {
  if (c.q && !(f.id.toLowerCase().includes(c.q) || (f.name || "").toLowerCase().includes(c.q)
      || (f.sublibrary || "").toLowerCase().includes(c.q))) return false;
  if (f.length != null && (f.length < c.lmin || f.length > c.lmax)) return false;
  if ((f.plddt || 0) < c.plddt) return false;
  if (f.clashscore != null && f.clashscore > c.clash) return false;
  if (c.tmhas && f.best_tm1 == null) return false;
  if (c.novelOnly && f.is_novel_v341 !== 1) return false;
  if (f.best_tm1 != null && f.best_tm1 > c.tmax) return false;
  if (f.overlap_ae != null && f.overlap_ae > c.ovmax) return false;
  if (c.shape && !f.shape_ok) return false;
  if (c.agr > -1 && (f.shape_agr == null || f.shape_agr < c.agr)) return false;
  if (c.crmin > 0 && (f.contact_ratio == null || f.contact_ratio < c.crmin)) return false;
  if (c.bpmin > 0 && (f.bp_fraction == null || f.bp_fraction < c.bpmin)) return false;
  if (c.tert && f.n_tert < 1) return false;
  if (c.rare && f.n_rare < 1) return false;
  if (c.motifs.length && !c.motifs.some((m) => (f.motifs || []).includes(m))) return false;
  if (c.pk !== "any" && String(f.pseudoknot) !== c.pk) return false;
  if (DSBYID[f._dsid] && DSBYID[f._dsid].motifs && f.letter && !c.letters.has(f.letter)) return false;
  return true;
}

function ranker(c) {
  const [key, dir] = (sortOverride ? `${sortOverride.key}:${sortOverride.dir}` : c.rank).split(":");
  const sign = dir === "asc" ? 1 : -1;
  return (a, b) => {
    let x = a[key], y = b[key];
    const xn = (x == null), yn = (y == null);
    if (xn && yn) return 0;
    if (xn) return 1;       // nulls always last
    if (yn) return -1;
    if (typeof x === "string") return sign * x.localeCompare(y);
    return sign * (x - y);
  };
}

function render() {
  const c = filters();
  let rows = FOLDS.filter((f) => pass(f, c));
  rows.sort(ranker(c));
  if (c.perLetter) {
    const seen = {};
    rows = rows.filter((f) => { seen[f.letter] = (seen[f.letter] || 0) + 1; return seen[f.letter] <= c.topn; });
  } else {
    rows = rows.slice(0, c.topn);
  }
  $("count").textContent = `${rows.length} shown`;
  drawTable(rows);
}

const COLS = [
  ["id", "seq_id"], ["name", "name"], ["letter", "L"], ["length", "len"], ["plddt", "pLDDT"],
  ["best_tm1", "best_tm1"], ["near", "nearest"], ["overlap_ae", "ovlp_AE"],
  ["shape_ok", "SHAPE"], ["shape_agr", "SHAPE agr"], ["contact_ratio", "compact"], ["bp_fraction", "paired"],
  ["n_tert", "tert"], ["n_rare", "rare"], ["pseudoknot", "PK"], ["motifs", "motifs"],
];

function drawTable(rows) {
  $("thead").innerHTML = COLS.map(([k, lbl]) => {
    let cls = "";
    if (sortOverride && sortOverride.key === k) cls = sortOverride.dir === "asc" ? "asc" : "sorted";
    return `<th data-k="${k}" class="${cls}">${lbl}</th>`;
  }).join("");
  $("thead").querySelectorAll("th").forEach((th) => th.addEventListener("click", () => {
    const k = th.dataset.k;
    if (sortOverride && sortOverride.key === k) sortOverride.dir = sortOverride.dir === "asc" ? "desc" : "asc";
    else sortOverride = { key: k, dir: "asc" };
    render();
  }));
  const body = rows.map((f) => {
    const chips = (f.motifs || []).slice(0, 6).map((m) =>
      `<span class="motif-chip" style="background:${motifColor(m)}">${m.replace(/_/g, " ").toLowerCase()}</span>`).join("");
    const pl = f.plddt || 0;
    const plbar = `<span class="bar" style="width:${pl * 0.45}px;background:${pl > 80 ? "#3a7d44" : pl > 60 ? "#edae49" : "#c0504d"}"></span> ${num(pl, 0)}`;
    const hasShape = f.r2a3 != null || f.mean_prot_2a3 != null;
    const shape = f.shape_ok ? `<span class="pill" style="background:#3a7d44">yes</span>`
      : (hasShape ? `<span class="muted" title="has SHAPE data but motif residues not protected">no</span>`
                  : `<span class="muted" title="no usable SHAPE data for this fold">n/d</span>`);
    return `<tr data-id="${esc(f.id)}">
      <td>${esc(f.id)}</td><td>${esc(f.name || "")}</td><td>${f.letter}</td>
      <td class="num">${f.length ?? ""}</td><td class="num">${plbar}</td>
      <td class="num">${num(f.best_tm1, 3)}</td><td title="${(f.near_title || "").replace(/"/g, "&quot;")}">${f.near || ""}</td>
      <td class="num">${num(f.overlap_ae, 2)}</td><td>${shape}</td>
      <td class="num">${num(f.shape_agr, 2)}</td><td class="num">${num(f.contact_ratio, 2)}</td><td class="num">${num(f.bp_fraction, 2)}</td>
      <td class="num">${f.n_tert}</td><td class="num">${f.n_rare}</td>
      <td>${f.pseudoknot ? "&#10003;" : ""}</td><td>${chips}</td></tr>`;
  }).join("");
  $("tbody").innerHTML = body;
  $("tbody").querySelectorAll("tr").forEach((tr) =>
    tr.addEventListener("click", () => openDeep(tr.dataset.id)));
}

// ---------------- deep view ----------------
let FBYK = {};
function foldById(id) { if (!FBYK[id]) FOLDS.forEach((f) => FBYK[f.id] = f); return FBYK[id]; }

function currentMode() { return localStorage.getItem("atlas_deepmode") || "modal"; }
function updateLayout() {
  const open = !$("deep").classList.contains("hidden");
  const m = currentMode();
  document.body.classList.remove("deepopen-right", "deepopen-bottom");
  if (open && m === "right") document.body.classList.add("deepopen-right");
  if (open && m === "bottom") document.body.classList.add("deepopen-bottom");
  document.querySelectorAll('#layoutctl button[data-mode]').forEach((b) => b.classList.toggle("active", b.dataset.mode === m));
}
function setDeepMode(m) {
  const d = $("deep");
  d.classList.remove("mode-modal", "mode-right", "mode-bottom");
  d.classList.add("mode-" + m);
  localStorage.setItem("atlas_deepmode", m);
  updateLayout();
  if (viewer) { try { viewer.resize(); viewer.render(); } catch (e) {} }
}
function closeDeep() {
  $("deep").classList.add("hidden");
  updateLayout();
  if (viewer) { try { viewer.resize(); } catch (e) {} }
}

async function openDeep(id) {
  const f = foldById(id);
  const ds = dsFor(f);
  $("deep").classList.remove("hidden");
  setDeepMode(currentMode());
  $("deep-title").textContent = `${id}${f.name ? "  —  " + f.name : ""}`;
  drawProps(f);
  let react = null;
  if (ds.react) {
    try { react = await (await fetch(durl(prefix(ds) + "react/" + id + ".json"))).json(); } catch (e) { react = null; }
  }
  currentDeep = { f, react };
  drawTracks(f, react);
  load3D(f, react);
}
let currentDeep = null;

function drawProps(f) {
  const verdict = f.best_tm1 == null ? "n/a (not scored vs v341)"
    : f.best_tm1 < 0.40 ? "novel" : f.best_tm1 < 0.45 ? "borderline" : "matches known fold";
  const rowsHtml = [
    ["Source", `${f.source} (lib ${f.letter})`],
    ["Sublibrary", f.sublibrary],
    ["Length", `${f.length} nt`],
    ["pLDDT / gpde", `${num(f.plddt)} / ${num(f.gpde, 3)}`],
    ["Clashscore", num(f.clashscore, 2)],
    ["Novelty (best_tm1 vs v341)", f.best_tm1 == null ? "&mdash;" : `${num(f.best_tm1, 3)} &mdash; ${verdict}`],
    ["Nearest known fold", `${f.near || "&mdash;"}${f.near_title ? ` &mdash; ${f.near_title}` : ""}`],
    ["Distinct vs A&ndash;E (overlap)", num(f.overlap_ae, 3)],
    ["SHAPE-supported", `${f.shape_ok ? "yes" : "no"} (SHAPE–pairing agreement = ${num(f.shape_agr, 3)}, + = good; mean prot = ${num(f.mean_prot_2a3, 3)})`],
    ["OpenKnot score", num(f.openknot, 3)],
    ["Pseudoknot", f.pseudoknot ? "yes" : "no"],
    ["Secondary-structure class", f.ss_class],
    ["Compactness (C1′ contact ratio)", num(f.contact_ratio, 3)],
    ["Base-paired fraction", num(f.bp_fraction, 3)],
    ["Tertiary motifs", `${f.n_tert} (rare ${f.n_rare})`],
  ];
  const chips = (f.motifs || []).map((m) =>
    `<span class="motif-chip" style="background:${motifColor(m)}">${m.replace(/_/g, " ").toLowerCase()}</span>`).join(" ");
  $("props").innerHTML = "<table>" + rowsHtml.map(([k, v]) => `<tr><td class="muted">${k}</td><td>${v}</td></tr>`).join("")
    + `<tr><td class="muted">Motifs</td><td>${chips}</td></tr></table>`;
}

function spansFor(f) {
  const M = MOTIFS_BY_DS[f._dsid] || {};
  return (M[f.id] || []).map(([type, res]) => {
    const ranges = [];
    res.split(",").forEach((c) => {
      const rng = c.trim().split(":").pop();
      const [a, b] = rng.includes("-") ? rng.split("-") : [rng, rng];
      const ai = parseInt(a), bi = parseInt(b);
      if (!Number.isNaN(ai)) ranges.push([ai, Number.isNaN(bi) ? ai : bi]);
    });
    return { type, ranges };
  });
}

function drawTracks(f, react) {
  const seq = (react && react.seq) || "";
  const ra = react && (react.a23 || react.dms);
  const n = seq.length || (ra ? ra.length : 0) || f.length || 0;
  if (!n) { $("tracks").innerHTML = '<p class="muted">No reactivity / sequence available for this fold.</p>'; return; }
  const cw = Math.max(6, Math.min(16, Math.floor(900 / n)));
  const W = n * cw, pad = 4;
  const motifs = spansFor(f);
  // lane assignment
  const lanes = [];
  motifs.forEach((m) => {
    const lo = Math.min(...m.ranges.map((r) => r[0])), hi = Math.max(...m.ranges.map((r) => r[1]));
    let li = lanes.findIndex((end) => end < lo);
    if (li < 0) { li = lanes.length; lanes.push(0); }
    lanes[li] = hi; m.lane = li;
  });
  const laneH = 9, mh = lanes.length * (laneH + 2);
  const yMot = pad, ySeq = yMot + mh + 4, yDms = ySeq + 16, yA23 = yDms + 16, yPair = yA23 + 16, H = yPair + 18;
  let svg = `<svg width="${W + 40}" height="${H}" font-size="9">`;
  // motif bars
  motifs.forEach((m) => {
    m.ranges.forEach(([a, b]) => {
      svg += `<rect x="${(a - 1) * cw}" y="${yMot + m.lane * (laneH + 2)}" width="${(b - a + 1) * cw}" height="${laneH}" rx="2" fill="${motifColor(m.type)}"><title>${m.type} ${a}-${b}</title></rect>`;
    });
  });
  // sequence
  for (let i = 0; i < n; i++) {
    const ch = seq[i] || "N";
    svg += `<rect x="${i * cw}" y="${ySeq}" width="${cw - 0.5}" height="13" fill="${nucColor(ch, altPalette())}"/>`;
    if (cw >= 10) svg += `<text x="${i * cw + cw / 2}" y="${ySeq + 10}" text-anchor="middle" fill="#fff">${ch}</text>`;
  }
  // reactivity rows
  const rrow = (arr, y, label) => {
    let s = `<text x="${W + 3}" y="${y + 11}" fill="#5b6670">${label}</text>`;
    for (let i = 0; i < n; i++) {
      const v = arr ? arr[i] : null;
      s += `<rect x="${i * cw}" y="${y}" width="${cw - 0.5}" height="13" fill="${shapeColor(v)}"><title>${label} ${i + 1}: ${v == null ? "n/a" : (+v).toFixed(2)}</title></rect>`;
    }
    return s;
  };
  svg += rrow(react && react.dms, yDms, "DMS");
  svg += rrow(react && react.a23, yA23, "2A3");
  // predicted pairing track: unpaired = light red, paired = white (eyeball SHAPE agreement)
  const dbn = (PAIRING_BY_DS[f._dsid] || {})[f.id] || "";
  let pr = `<text x="${W + 3}" y="${yPair + 11}" fill="#5b6670">pair</text>`;
  for (let i = 0; i < n; i++) {
    const ch = dbn[i];
    const paired = ch && ch !== "." && ch !== "-";
    const fill = ch ? (paired ? "#ffffff" : "#f3a0a0") : "#eef2f5";
    pr += `<rect x="${i * cw}" y="${yPair}" width="${cw - 0.5}" height="13" fill="${fill}" stroke="#dfe3e8" stroke-width="0.5"><title>pos ${i + 1}: ${ch ? (paired ? "paired" : "unpaired") : "n/a"}</title></rect>`;
  }
  svg += pr + "</svg>";
  $("tracks").innerHTML = `<div style="font-size:11px;color:#5b6670;margin-bottom:3px">motif lanes &middot; sequence &middot; DMS &middot; 2A3 reactivity (white=protected &rarr; red=reactive) &middot; pairing (white=paired, light red=unpaired)</div>` + svg;
}

let viewer = null;
async function load3D(f, react) {
  const el = $("viewer3d");
  el.innerHTML = "";
  if (typeof $3Dmol === "undefined") { el.innerHTML = '<p style="color:#fff;padding:8px">3Dmol.js not loaded.</p>'; return; }
  const id = f.id;
  const ds = dsFor(f);
  let data;
  try { data = await (await fetch(durl(prefix(ds) + "structs/" + (f.key || f.id) + "." + (ds.ext || "cif")))).text(); }
  catch (e) { el.innerHTML = '<p style="color:#fff;padding:8px">structure unavailable</p>'; return; }
  viewer = $3Dmol.createViewer(el, { backgroundColor: "0x0d1117" });
  const fmt = data.startsWith("data_") || data.includes("_atom_site") ? "cif" : "pdb";
  viewer.addModel(data, fmt);
  const a23 = react && react.a23, dms = react && react.dms;
  const dbn = (PAIRING_BY_DS[f._dsid] || {})[f.id] || "";
  const mode = ($("color_by") && $("color_by").value) || "a23";
  let style;
  if (mode === "plddt") {                       // per-residue pLDDT from the B-factor (AlphaFold palette)
    style = { cartoon: { colorfunc: (a) => plddtCF(a.b), ringMode: 3 } };
  } else if (mode === "pairing") {              // predicted secondary structure: paired vs unpaired
    style = { cartoon: { colorfunc: (a) => { const c = dbn[a.resi - 1]; const p = c && c !== "." && c !== "-"; return p ? "0xffffff" : (c ? "0xf3a0a0" : "0x9aa7b0"); }, ringMode: 3 } };
  } else if (mode === "nuc") {                   // nucleotide identity (same palette as the sequence track)
    const alt = altPalette();
    style = { cartoon: { colorfunc: (a) => hexCF(nucColor((a.resn || "").trim(), alt)), ringMode: 3 } };
  } else if (mode === "spectrum") {
    style = { cartoon: { color: "spectrum", ringMode: 3 } };
  } else {                                       // a23 / dms reactivity (blue protected -> red reactive); spectrum if absent
    const arr = mode === "dms" ? dms : a23;
    style = arr ? { cartoon: { colorfunc: (a) => reactCF(arr[a.resi - 1]), ringMode: 3 } }
                : { cartoon: { color: "spectrum", ringMode: 3 } };
  }
  viewer.setStyle({}, style);
  spansFor(f).forEach((m) => {
    const resi = [];
    m.ranges.forEach(([a, b]) => { for (let r = a; r <= b; r++) resi.push(r); });
    viewer.addStyle({ resi }, { stick: { color: motifColor(m.type), radius: 0.28 } });
  });
  viewer.zoomTo(); viewer.render();
}
function shapeMix(a, b, t) {
  return "0x" + [0, 1, 2].map((i) => Math.round((a[i] + (b[i] - a[i]) * t) * 255).toString(16).padStart(2, "0")).join("");
}
function hexCF(h) { return "0x" + h.replace("#", ""); }
function reactCF(v) {                              // 2A3/DMS: blue protected -> white -> red reactive
  if (v == null || Number.isNaN(v)) return "0xf7f7f7";
  const t = Math.max(-0.3, Math.min(1, v)), f = (t + 0.3) / 1.3;
  return f < 0.5 ? $3Dmol.CC.color(shapeMix(PROT_STOPS[0], PROT_STOPS[1], f * 2))
                 : $3Dmol.CC.color(shapeMix(PROT_STOPS[1], PROT_STOPS[2], (f - 0.5) * 2));
}
function plddtCF(b) {                              // AlphaFold confidence palette
  if (b == null || Number.isNaN(b)) return "0xcccccc";
  return b >= 90 ? "0x0053d6" : b >= 70 ? "0x65cbf3" : b >= 50 ? "0xffdb13" : "0xff7d45";
}

boot();
