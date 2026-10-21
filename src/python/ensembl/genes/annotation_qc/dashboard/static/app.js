"use strict";
// Annotation QC dashboard (no external dependencies). Every number shown is a stored
// comparator metric row (or a legacy-import row); this script selects and formats them.
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const esc = (v) => String(v ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmtInt = (n) => (n === null || n === undefined || Number.isNaN(n) ? "—" : Math.round(n).toLocaleString("en-GB"));
const fmtPct = (v, d = 1) => (v === null || v === undefined ? "—" : (100 * v).toFixed(d) + "%");
const ALL = "all";
const DISPLAYED = ["complete", "unverified", "legacy"];
const RESULT_LABEL = { complete: "complete", unverified: "unverified implementation", legacy: "pre-contract (legacy import)", failed: "failed",
  not_run: "not run", stale: "stale", unavailable: "unavailable" };
const OUTCOME = {
  exact: ["exact", "coordinate-exact CDS (best CDS pair)"],
  terminal_diff: ["ends ≠", "same CDS intron chain; start/stop differs"],
  diff_chain: ["splice ≠", "≥ 0.8 reciprocal CDS overlap, different CDS intron chain"],
  partial_lt08: ["< 0.8", "shares CDS, < 0.8 reciprocal overlap"],
  span_only: ["span", "same-strand gene-span partner without CDS overlap"],
  strand_mismatch: ["strand", "opposite-strand overlap only"],
  missed: ["missed", "no partner"],
};
const HEADLINE = [
  { id: "locus_recovery_cds_overlap", label: "Loci recovered", sub: "same-strand CDS overlap · % of R", fmt: "pct" },
  { id: "exact_cds_one_to_one_recall", label: "1:1 exact CDS recall", sub: "TP / R", fmt: "pct" },
  { id: "exact_cds_one_to_one_precision", label: "1:1 exact CDS precision", sub: "TP / Q", fmt: "pct" },
  { id: "exact_cds_one_to_one_f1", label: "1:1 exact CDS F1", sub: "2TP / (R + Q) · 0–1", fmt: "f1" },
];
const STRUCTURE = [
  { id: "cds_intron_chain_best_pair", label: "CDS intron chain identical", sub: "best pair · % of multi-segment R" },
  { id: "cds_intron_chain_any_pair", label: "CDS intron chain identical", sub: "any pair · % of multi-segment R" },
  { id: "ref_unique_cds_introns_recovered", label: "Reference CDS introns recovered", sub: "unique introns" },
  { id: "query_unique_cds_introns_supported", label: "Predicted CDS introns supported", sub: "unique introns" },
  { id: "split_reference_genes", label: "Reference genes split", sub: "% of R" },
  { id: "merge_query_genes", label: "Predictions merging genes", sub: "% of Q" },
];
const NOVEL = ["query_novel_Overlaps_reference_pseudogene", "query_novel_Overlaps_other_non_protein_coding_reference",
  "query_novel_Overlaps_unevaluated_protein_coding_reference", "query_novel_Overlaps_evaluated_reference_opposite_strand",
  "query_novel_Novel_no_reference_overlap"];
const S = { exp: null, meta: null, ref: null, policy: null, scope: ALL, view: "headline", showSupp: true, ov: null, idx: new Map(), dict: {},
  genePage: 0, predPage: 0, pageSize: 100, locus: null, selectedGene: null, posCoord: "cds_span", posMin: 0, lenMeasure: "exact", splice: {} };

async function api(path, params = {}) {
  const res = await fetch(apiUrl(path, params));
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || res.statusText);
  return data;
}
function apiUrl(path, params = {}) {
  const url = new URL(path, location.origin);
  for (const [k, v] of Object.entries({ exp: S.exp, ...params })) if (v !== undefined && v !== null && v !== "") url.searchParams.set(k, v);
  return url.toString();
}
function showError(where, err) { where.innerHTML = `<div class="issue error">${esc(err.message || err)}</div>`; }

// ---------------------------------------------------------------- tooltip (mouse and keyboard focus)
const tip = $("#tooltip");
function showTip(el, x, y) { tip.textContent = el.dataset.tip; tip.hidden = false; tip.style.left = Math.min(x + 14, window.innerWidth - 390) + "px"; tip.style.top = (y + 14) + "px"; }
document.addEventListener("mousemove", (e) => { const t = e.target.closest("[data-tip]"); if (!t) { tip.hidden = true; return; } showTip(t, e.clientX, e.clientY); });
document.addEventListener("focusin", (e) => { const t = e.target.closest("[data-tip]"); if (!t) { tip.hidden = true; return; } const r = t.getBoundingClientRect(); showTip(t, r.left, r.bottom); });
document.addEventListener("keydown", (e) => { if (e.key === "Escape") tip.hidden = true; });

// ---------------------------------------------------------------- keyboard: tab lists and radio groups
function arrowNav(container, selector, activate) {
  container.addEventListener("keydown", (e) => {
    if (!["ArrowLeft", "ArrowRight", "Home", "End"].includes(e.key)) return;
    const items = $$(selector, container); const i = items.indexOf(document.activeElement);
    if (i < 0) return;
    e.preventDefault();
    const j = e.key === "Home" ? 0 : e.key === "End" ? items.length - 1 : (i + (e.key === "ArrowRight" ? 1 : -1) + items.length) % items.length;
    items[j].focus(); activate(items[j]);
  });
}
function setTabs(container, attr, value) {
  $$("button", container).forEach((b) => { const on = b.dataset[attr] === value; b.setAttribute("aria-selected", on); b.tabIndex = on ? 0 : -1; });
}

// ---------------------------------------------------------------- setup
async function init() {
  const exps = await api("/api/experiments");
  const species = (e) => (e.species ? e.species.charAt(0).toUpperCase() + e.species.slice(1) : e.title || e.id);
  $("#exp").innerHTML = exps.map((e) => `<option value="${esc(e.id)}" title="${esc(e.title || e.id)}">${esc(species(e))}${e.synthetic ? " [SYNTHETIC]" : ""}</option>`).join("");
  const params = new URLSearchParams(location.hash.slice(1));
  S.exp = exps.some((e) => e.id === params.get("exp")) ? params.get("exp") : exps[0].id;
  $("#exp").value = S.exp;
  $("#exp").onchange = () => { S.exp = $("#exp").value; loadExperiment(); };
  $$(".tabs button").forEach((b) => (b.onclick = () => selectTab(b.dataset.tab)));
  $$(".subnav button").forEach((b) => (b.onclick = () => selectView(b.dataset.view)));
  arrowNav($(".tabs"), "button", (b) => selectTab(b.dataset.tab));
  arrowNav($(".subnav"), "button", (b) => selectView(b.dataset.view));
  $("#showSupp").onchange = () => { S.showSupp = $("#showSupp").checked; renderAll(); };
  await loadExperiment(params);
}
function selectView(name) {
  S.view = name; setTabs($(".subnav"), "view", name);
  $$(".view").forEach((v) => v.classList.toggle("active", v.id === "view-" + name));
  renderView(); writeHash();
}
function selectTab(name) {
  setTabs($(".tabs"), "tab", name);
  $$(".tab").forEach((t) => t.classList.toggle("active", t.id === "tab-" + name));
  if (name === "genes" && !S.genesLoaded) loadGenes();
  if (name === "predictions" && !S.predsLoaded) loadPredictions();
  writeHash();
}
async function loadExperiment(params = new URLSearchParams()) {
  S.meta = await api("/api/meta");
  S.dict = Object.fromEntries((S.meta.metric_dictionary?.metrics || []).map((m) => [m.id, m]));
  $("#synthetic").hidden = !S.meta.experiment.synthetic;
  const c = S.meta.contract || {};
  $("#contractNote").textContent = `comparator output contract ${c.contract_version || "?"} · definitions ${c.definition_version || "?"}`;
  $("#footContract").textContent = `${c.contract_id || ""} ${c.contract_version || ""}`;
  // One species reference per experiment (`annotation-qc dashboard` refuses other configurations).
  S.ref = S.meta.references[0].id;
  const pols = S.meta.policies.map((p) => p.id);
  S.policy = pols.includes(params.get("policy")) ? params.get("policy") : pols.includes("A") ? "A" : pols[0];
  S.scope = params.get("scope") || ALL;
  if (["cds_span", "coding"].includes(params.get("pos"))) S.posCoord = params.get("pos");
  if (["0", "60"].includes(params.get("posmin"))) S.posMin = Number(params.get("posmin"));
  if (["exact", "donor", "acceptor"].includes(params.get("len"))) S.lenMeasure = params.get("len");
  if (params.get("supp") === "0") { S.showSupp = false; $("#showSupp").checked = false; }
  S.splice = {};
  renderDictionary();
  await loadReference();
  if (["headline", "structure", "detail"].includes(params.get("view"))) selectView(params.get("view"));
  if (["genes", "predictions", "provenance", "dictionary"].includes(params.get("tab"))) selectTab(params.get("tab"));
}
function refMeta() { return S.meta.references.find((r) => r.id === S.ref); }
function scopes() { return refMeta().scopes || [{ id: ALL, label: "All compared sequences", available: true }]; }
function scopeInfo() { return scopes().find((x) => x.id === S.scope) || scopes()[0]; }
function policyInfo(id = S.policy) { return S.meta.policies.find((p) => p.id === id); }
function runById(id) { return S.meta.runs.find((r) => r.run_id === id); }
function configRuns() { return S.meta.runs.filter((r) => r.reference === S.ref); }
function result(run, policy = S.policy) { return S.ov.results.find((r) => r.run_id === run && r.policy_id === policy); }
function shown(res) { return res && DISPLAYED.includes(res.status); }
function metric(run, id, policy = S.policy, scope = S.scope) { return S.idx.get(`${run}|${policy}|${scope}|${id}`); }
function runName(r) { return r.label || r.tool; }
// Display roles (configuration runs[].display_role): presentation only, never a value.
const ROLE_STYLE = {
  focus: { rank: 0, tag: "focus", width: 3.5, r: 5, shape: "circle" },
  primary_comparator: { rank: 1, tag: "primary comparison", width: 2.75, r: 4.5, shape: "square" },
};
const OTHER_STYLE = { rank: 2, tag: "", width: 1.5, r: 3.5, shape: "ring" };
function roleStyle(r) { return ROLE_STYLE[r.display_role] || OTHER_STYLE; }
function marker(shape, x, y, r, color, tip) {
  const t = tip ? ` data-tip="${esc(tip)}"` : "";
  if (shape === "square") return `<rect x="${(x - r).toFixed(1)}" y="${(y - r).toFixed(1)}" width="${2 * r}" height="${2 * r}" fill="${color}" stroke="#fff" stroke-width="1"${t}/>`;
  if (shape === "ring") return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${r}" fill="#fff" stroke="${color}" stroke-width="1.75"${t}/>`;
  return `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${r}" fill="${color}" stroke="#fff" stroke-width="1.25"${t}/>`;
}
function swatch(r) { const st = roleStyle(r); return `<svg class="mark" width="14" height="14" viewBox="0 0 14 14" aria-hidden="true">${marker(st.shape, 7, 7, st.shape === "ring" ? 4.5 : 5.5, r.color)}</svg>`; }
function roleTag(r) { const st = ROLE_STYLE[r.display_role]; return st ? `<span class="tag role role-${esc(r.display_role)}">${esc(st.tag)}</span>` : ""; }
function isSupp(r) { return r.role === "supplementary" || !!r.companion_of; }
function harmonised(r) { return r.variant === "stop_harmonised"; }

// Order: runs with a display role (focus, then primary comparator), pinned runs, preferred tools (configured
// order), then the other primary runs by descending display.order_metric (default one-to-one exact-CDS recall)
// in the selected policy and scope. A supplementary companion follows its run and is never ranked.
function orderRuns(runs) {
  const d = S.meta.display || {};
  const byConfig = (a, b) => a.config_order - b.config_order;
  const value = (r) => { const m = metric(r.run_id, d.order_metric || "exact_cds_one_to_one_recall"); return m && m.status === "ok" && shown(result(r.run_id)) ? m.value : -1; };
  const companion = (r) => r.companion_of && runs.some((x) => x.run_id === r.companion_of);
  const leaders = runs.filter((r) => !companion(r));
  const roled = leaders.filter((r) => ROLE_STYLE[r.display_role] && !isSupp(r))
    .sort((a, b) => roleStyle(a).rank - roleStyle(b).rank || value(b) - value(a) || byConfig(a, b));
  const pinned = (d.pinned_runs || []).map((id) => leaders.find((r) => r.run_id === id)).filter((r) => r && !roled.includes(r));
  const preferred = [];
  for (const tool of d.preferred_tools || []) leaders.filter((r) => r.tool === tool && !roled.includes(r) && !pinned.includes(r)).sort((a, b) => value(b) - value(a) || byConfig(a, b)).forEach((r) => preferred.push(r));
  const rest = leaders.filter((r) => !roled.includes(r) && !pinned.includes(r) && !preferred.includes(r))
    .sort((a, b) => (isSupp(a) - isSupp(b)) || value(b) - value(a) || byConfig(a, b));
  const out = [];
  for (const r of [...roled, ...pinned, ...preferred, ...rest]) {
    out.push(r);
    runs.filter((c) => companion(c) && c.companion_of === r.run_id).sort(byConfig).forEach((c) => out.push(c));
  }
  return out;
}
function visibleRuns() { return orderRuns(configRuns()).filter((r) => S.showSupp || !isSupp(r)); }

async function loadReference() {
  S.ov = await api("/api/overview", { reference: S.ref });
  S.idx = new Map(S.ov.metrics.map((m) => [`${m.run_id}|${m.policy_id}|${m.scope_id}|${m.metric_id}`, m]));
  if (!scopes().some((x) => x.id === S.scope && x.available)) S.scope = ALL;
  S.genesLoaded = S.predsLoaded = false; S.splice = {};
  $("#locusPanel").hidden = true;
  setupGeneFilters(); setupPredFilters(); renderAll(); renderProvenance();
  const active = $(".tabs button[aria-selected=true]").dataset.tab;
  if (active === "genes") loadGenes();
  if (active === "predictions") loadPredictions();
}
function renderAll() { $("#showSupp").closest("label").hidden = !configRuns().some(isSupp); renderPolicySelector(); renderScopeSelector(); renderContextBar(); renderView(); writeHash(); }
function renderView() {
  if (S.view === "headline") renderHeadline();
  if (S.view === "structure") { renderStructure(); renderPlots(); }
  if (S.view === "detail") renderDetail();
}
function writeHash() {
  const tab = $(".tabs button[aria-selected=true]")?.dataset.tab || "overview";
  history.replaceState(null, "", `#exp=${encodeURIComponent(S.exp)}&policy=${encodeURIComponent(S.policy)}`
    + `&scope=${encodeURIComponent(S.scope)}&view=${S.view}&tab=${tab}${S.showSupp ? "" : "&supp=0"}`
    + (S.view === "structure" ? `&pos=${S.posCoord}&posmin=${S.posMin}&len=${S.lenMeasure}` : ""));
}
function policyAvailability(pid) {
  const rs = configRuns().map((r) => result(r.run_id, pid)).filter(Boolean);
  if (!rs.length) return { ok: false, why: "no results" };
  if (!rs.some(shown)) return { ok: false, why: rs[0].reason || rs[0].status };
  return { ok: true };
}
function renderPolicySelector() {
  const box = $("#policy");
  box.innerHTML = S.meta.policies.map((p) => {
    const av = policyAvailability(p.id), on = p.id === S.policy;
    return `<button type="button" role="radio" aria-checked="${on}" tabindex="${on ? 0 : -1}" data-p="${esc(p.id)}"
      data-tip="${esc(`${p.label}\n${p.description || ""}${av.ok ? "" : "\nNo results: " + av.why}`)}">${esc(p.id)} · ${esc(p.label)}${av.ok ? "" : " (none)"}</button>`;
  }).join("");
  const pick = (b) => { S.policy = b.dataset.p; S.genesLoaded = S.predsLoaded = false; renderAll();
    const t = $(".tabs button[aria-selected=true]").dataset.tab; if (t === "genes") loadGenes(); if (t === "predictions") loadPredictions(); if (S.locus) drawLocus(); };
  $$("button", box).forEach((b) => (b.onclick = () => pick(b)));
  if (!box.dataset.nav) { arrowNav(box, "button", pick); box.dataset.nav = "1"; }
}
function renderScopeSelector() {
  const sel = $("#scope");
  sel.innerHTML = scopes().map((x) => `<option value="${esc(x.id)}" ${x.available ? "" : "disabled"}>${esc(x.label)}${x.available ? "" : " — unavailable"}</option>`).join("");
  sel.value = S.scope; sel.disabled = scopes().length < 2;
  sel.onchange = () => { S.scope = sel.value; S.genesLoaded = S.predsLoaded = false; renderAll();
    const t = $(".tabs button[aria-selected=true]").dataset.tab; if (t === "genes") loadGenes(); if (t === "predictions") loadPredictions(); };
}
function identitySummary() {
  const counts = {};
  for (const r of S.ov.results) if (shown(r)) {
    const id = r.identity || {};
    const key = r.status === "legacy" ? "pre-contract results (no implementation fingerprint)" : id.status === "verified" ? `verified implementation ${id.build_id} (historical build)`
      : id.status === "current_checkout" ? "current code (this checkout)" : "unverified implementation";
    counts[key] = (counts[key] || 0) + 1;
  }
  return Object.entries(counts).map(([k, n]) => `${n} × ${k}`).join("; ") || "no results";
}
function definitionSummary() {
  const v = [...new Set(S.ov.results.filter(shown).map((r) => r.definition_version || (r.status === "legacy" ? "pre-contract" : "?")))];
  return v.join(", ") || "—";
}
function renderContextBar() {
  const p = policyInfo(), sc = scopeInfo(), ref = refMeta();
  const R = configRuns().map((r) => metric(r.run_id, "reference_genes")).find((m) => m && m.status === "ok");
  const supp = configRuns().filter(isSupp);
  $("#contextBar").innerHTML = `<span class="ctx"><b>Experiment</b>${esc(S.meta.experiment.title)}</span>
    <span class="ctx"><b>Reference</b>${esc(ref.species)} · ${esc(ref.assembly)}${ref.assembly_accession ? " (" + esc(ref.assembly_accession) + ")" : ""} · ${esc(ref.annotation_release || "")}</span>
    <span class="ctx" data-tip="${esc(p.description || "")}" tabindex="0"><b>Transcript policy</b>${esc(p.id)} · ${esc(p.label)}</span>
    <span class="ctx" data-tip="${esc(sc.description || "")}" tabindex="0"><b>Sequence scope</b>${esc(sc.label)}${sc.sequence_count ? ` (${fmtInt(sc.sequence_count)} sequences)` : ""} · R = ${R ? fmtInt(R.numerator) : "—"}</span>
    ${supp.length ? `<span class="ctx"><b>Supplementary variants</b>${supp.length} ${S.showSupp ? "shown beside their run" : "hidden"}</span>` : ""}
    <span class="ctx" data-tip="${esc(`Comparator output contract ${S.meta.contract?.contract_version || "?"}; definitions of the displayed results: ${definitionSummary()}`)}" tabindex="0"><b>Comparison</b>definitions ${esc(definitionSummary())}</span>
    <span class="ctx"><b>Implementation</b>${esc(identitySummary())}</span>`;
}

// ---------------------------------------------------------------- cells
const NA_TEXT = { undefined: "undefined", not_applicable: "not applicable", not_available: "not available" };
function naCell(text, why, cls = "num") { return `<td class="${cls}"><span class="na" tabindex="0" data-tip="${esc(why || text)}">${esc(text)}</span></td>`; }
function resultGate(run) {
  const res = result(run.run_id);
  if (!res) return { ok: false, text: "no result", why: "no result for this run and policy" };
  if (!shown(res)) return { ok: false, text: RESULT_LABEL[res.status] || res.status, why: res.reason };
  return { ok: true, res };
}
function fmtValue(m, fmt) {
  const v = m.value;
  if (fmt === "f1") return v > 0 && v < 0.001 ? "< 0.001" : v.toFixed(3);
  if (fmt === "count") return fmtInt(m.numerator);
  return v > 0 && v < 0.001 ? "< 0.1%" : fmtPct(v);
}
function metricTip(m, id) {
  const def = (S.meta.metric_definitions && Object.values(S.meta.metric_definitions)[0] || {})[id] || {};
  const parts = [`${S.dict[id]?.label || id}`, `numerator ${fmtInt(m.numerator)} (${m.numerator_id || def.numerator || ""})`];
  if (m.denominator !== null && m.denominator !== undefined) parts.push(`denominator ${fmtInt(m.denominator)} (${m.denominator_id || ""})`);
  if (m.value !== null && m.value !== undefined && m.denominator !== null && m.denominator !== undefined) parts.push(`value ${m.value_text || m.value}`);
  parts.push(`status ${m.status}${m.reason ? " — " + m.reason : ""} · source ${m.source}`);
  return parts.join("\n");
}
function valueCell(run, id, fmt, opts = {}) {
  const gate = resultGate(run);
  if (!gate.ok) return naCell(gate.text, gate.why);
  const m = metric(run.run_id, id);
  if (!m) return naCell("not available", "no row for this metric, policy and scope in the stored results");
  if (m.status !== "ok") return naCell(NA_TEXT[m.status] || m.status, `${m.status}: ${m.reason || "no reason recorded"}`);
  const v = m.value;
  const sub = opts.sub === false ? "" : `<div class="cell-sub">${m.denominator !== null && m.denominator !== undefined ? `${fmtInt(m.numerator)} / ${fmtInt(m.denominator)}` : ""}</div>`;
  const bar = opts.bar === false || v === null ? "" : `<div class="minibar"><span style="width:${(100 * Math.min(1, Math.max(0, v))).toFixed(1)}%;background:${run.color}"></span></div>`;
  const flag = gate.res.status === "unverified" ? ' <span class="tag bad" data-tip="unverified implementation">!</span>' : gate.res.status === "legacy" ? ' <span class="tag warn" data-tip="pre-contract result (legacy import)">pre</span>' : "";
  return `<td class="num" tabindex="0" data-tip="${esc(metricTip(m, id))}"><div class="cell-main">${fmtValue(m, fmt)}${flag}</div>${sub}${bar}</td>`;
}
function countCell(run, id, field = "numerator") {
  const gate = resultGate(run);
  if (!gate.ok) return naCell(gate.text, gate.why);
  const m = metric(run.run_id, id);
  if (!m) return naCell("not available", "no row for this metric, policy and scope");
  if (m[field] === null || m[field] === undefined) return naCell(NA_TEXT[m.status] || m.status, `${m.status}: ${m.reason || ""}`);
  return `<td class="num" tabindex="0" data-tip="${esc(metricTip(m, id))}">${fmtInt(m[field])}</td>`;
}
function runLimitations(r) {
  const res = result(r.run_id);
  const out = (res?.limitations || []).map((l) => l.text);
  for (const c of r.caveats || []) out.push(c);
  const ia = r.input_audit;
  if (ia && ia.available && ia.reference_pc_genes_on_sequences_without_records?.genes) {
    const pr = ia.processed_sequences_record || {};
    out.push(`Processing scope (input audit, whole assembly): predictions on ${fmtInt(ia.sequences_with_records)} sequences; ${pr.sequences ? `${fmtInt(pr.sequences.length ?? pr.sequences)} sequences recorded as processed (${pr.source})` : "which sequences were processed is not recorded"}. ${fmtInt(ia.reference_pc_genes_on_sequences_without_records.genes)} reference protein-coding genes on ${fmtInt(ia.reference_pc_genes_on_sequences_without_records.sequences)} sequences without prediction records count as missed; that does not establish that those sequences were left unprocessed.`);
  }
  return out;
}
function runCell(r, extra = "") {
  const res = result(r.run_id);
  const lim = runLimitations(r);
  const tags = [
    isSupp(r) ? `<span class="tag info" data-tip="${esc("Supplementary variant of " + (runById(r.companion_of) ? runName(runById(r.companion_of)) : "the original run") + ": not ranked as a separate tool")}">supplementary</span>` : "",
    roleTag(r),
    harmonised(r) && isSupp(r) ? `<span class="tag warn" data-tip="${esc("query compared after --add-stop-codon (genetic code " + (r.comparator_options?.stop_codon_genetic_code ?? "?") + (r.comparator_options?.stop_codon_sequence_code ? ", " + r.comparator_options.stop_codon_sequence_code.join(", ") : "") + "): terminal CDS extended by an adjacent in-frame stop where the frame is established from the sequence")}">stop-harmonised</span>` : "",
    res?.status === "unverified" ? `<span class="tag bad" data-tip="${esc(res.reason || "")}">unverified implementation</span>` : "",
    res?.status === "legacy" ? `<span class="tag warn" data-tip="${esc(res.reason || "")}">pre-contract</span>` : "",
  ].join(" ");
  const detail = [`run id ${r.run_id}`, r.tool !== runName(r) ? `tool ${r.tool}` : "", r.model ? `model ${r.model}` : "", r.tool_version ? `version ${r.tool_version}` : "",
    `annotation ${r.annotation.split("/").slice(-1)[0]}`, harmonised(r) ? `comparator option --add-stop-codon (genetic code ${r.comparator_options?.stop_codon_genetic_code ?? "?"}): terminal stop codon added where the file's CDS excludes it (audited)` : "",
    ROLE_STYLE[r.display_role] ? `display role: ${ROLE_STYLE[r.display_role].tag} (presentation only)` : ""].filter(Boolean).join("\n");
  const cav = lim.length ? ` <span class="caveat" tabindex="0" data-tip="${esc("Limitations:\n• " + lim.join("\n• "))}" aria-label="limitations">ⓘ</span>` : "";
  return `<td class="runcell${ROLE_STYLE[r.display_role] ? " role-" + esc(r.display_role) : ""}" style="--run:${r.color}">${isSupp(r) ? '<span class="muted" aria-hidden="true">↳ </span>' : ""}${swatch(r)}<span tabindex="0" data-tip="${esc(detail)}"><b>${esc(runName(r))}</b></span>${cav} ${tags}${extra}</td>`;
}

// ---------------------------------------------------------------- headline
function renderHeadline() {
  const runs = visibleRuns();
  const bad = runs.filter((r) => !resultGate(r).ok);
  const legacyRuns = runs.filter((r) => result(r.run_id)?.status === "legacy");
  $("#resultStatus").innerHTML = (bad.length ? `<div class="issue warning">${bad.map((r) => `${esc(runName(r))}: ${esc(resultGate(r).text)}${resultGate(r).why ? " — " + esc(resultGate(r).why) : ""}`).join("<br>")}</div>` : "")
    + (legacyRuns.length ? `<div class="issue info"><b>Pre-contract results (legacy import).</b> Only summary fields whose definition equals a contract metric are shown; metrics these results do not contain (e.g. the comparator's one-to-one matching or locus recovery, when they predate it) are “not available”, never estimated. The earlier best-pair estimates are listed separately under All metrics → Historical estimates and never fill the one-to-one columns.</div>` : "");
  const head = `<tr><th scope="col">Prediction run</th>${HEADLINE.map((c) => `<th scope="col" class="num"><span class="info" tabindex="0" data-tip="${esc(S.dict[c.id]?.definition || "")}">${esc(c.label)}</span><div class="cell-sub">${esc(c.sub)}</div></th>`).join("")}</tr>`;
  const tagged = (r, c) => valueCell(r, c.id, c.fmt).replace("<td ", `<td data-run="${esc(r.run_id)}" data-metric="${c.id}" `);
  const body = runs.map((r) => `<tr class="${isSupp(r) ? "supp" : ""}" data-run="${esc(r.run_id)}">${runCell(r)}${HEADLINE.map((c) => tagged(r, c)).join("")}</tr>`).join("");
  $("#summaryTable").innerHTML = runs.length ? `<div class="xscroll"><table class="grid summary"><caption class="sr">Headline CDS metrics, policy ${esc(S.policy)}, scope ${esc(scopeInfo().label)}</caption>${head}${body}</table></div>
    <p class="muted small">${esc(orderText(runs))}</p>` : '<p class="na">No prediction runs configured for this reference.</p>';
  renderOutcomes(runs);
  renderLimitBox(runs);
  const tpCols = [["reference_genes", "R reference genes"], ["query_genes", "Q predicted genes"], ["locus_recovery_cds_overlap", "loci recovered"], ["exact_cds_one_to_one_tp", "TP one-to-one pairs"],
    ["exact_cds_gene_pairs", "identical-CDS gene pairs (edges)"]];
  const pgCols = [["cds_coordinate_exact_reference", "reference genes coordinate-exact"], ["cds_coordinate_exact_query", "predicted genes coordinate-exact"],
    ["exact_cds_reference_with_candidate_unmatched", "reference genes with an identical-CDS partner, unmatched"], ["exact_cds_query_with_candidate_unmatched", "predicted genes with an identical-CDS partner, unmatched"]];
  $("#summaryCounts").innerHTML = `<div class="xscroll"><table class="grid compact"><tr><th rowspan="2">Prediction run</th><th colspan="${tpCols.length}" class="num">One-to-one exact CDS matching</th><th colspan="${pgCols.length}" class="num">Per-gene counts (best CDS pair or graph; not one-to-one)</th></tr>
    <tr>${tpCols.map(([, l]) => `<th class="num">${esc(l)}</th>`).join("")}${pgCols.map(([id, l]) => `<th class="num"><span class="info" tabindex="0" data-tip="${esc(S.dict[id]?.definition || "")}">${esc(l)}</span></th>`).join("")}</tr>
    ${runs.map((r) => `<tr class="${isSupp(r) ? "supp" : ""}">${runCell(r)}${tpCols.map(([id]) => countCell(r, id)).join("")}${pgCols.map(([id]) => countCell(r, id)).join("")}</tr>`).join("")}</table></div>
    <p class="muted small">Recall = TP / R; precision = TP / Q; F1 = 2TP / (R + Q). The per-gene coordinate-exact counts use each gene's best CDS pair among gene-span partners and can exceed TP when one gene has exact partners in several genes. ${S.scope === ALL ? "" : `Scope “${esc(scopeInfo().label)}”: the comparator's own rows for this scope (pairing, matching and introns never cross sequences).`}</p>`;
  const p = policyInfo();
  $("#policyExplain").textContent = p.reference_transcript_selection === "all"
    ? `Policy ${p.id}: a reference gene counts when any of its eligible isoforms is reproduced (any-isoform recovery).`
    : `Policy ${p.id}: one transcript per reference gene (${p.reference_transcript_selection.replace("_", " ")}) and per predicted gene (${p.query_transcript_selection.replace("_", " ")}); a correct model of another isoform counts as a miss.`;
}
function orderText(runs) {
  const d = S.meta.display || {}, metricLabel = S.dict[d.order_metric || "exact_cds_one_to_one_recall"]?.label || "one-to-one exact CDS recall";
  const first = [...runs.filter((r) => ROLE_STYLE[r.display_role] && !isSupp(r)).map((r) => `${runName(r)} (${ROLE_STYLE[r.display_role].tag})`),
    ...(d.preferred_tools || []).filter((t) => runs.some((r) => r.tool === t && !ROLE_STYLE[r.display_role]))];
  return `Order: ${first.length ? first.join(", ") + " first, then the other runs" : "runs"} by ${metricLabel} under this policy and scope (display order only).${runs.some(isSupp) ? " Supplementary variants follow their run and are not ranked." : ""}`;
}
// Reference genes by CDS outcome: the comparator's cds_outcome_* rows (they partition R); drawn, not recomputed.
const OUTCOME_ORDER = ["exact", "terminal_diff", "diff_chain", "partial_lt08", "span_only", "strand_mismatch", "missed"];
const OUTCOME_FILL = { exact: ["#1b7837", "#fff"], terminal_diff: ["#7fbc41", "#10290a"], diff_chain: ["#d9c45a", "#2b2405"], partial_lt08: ["#e08214", "#2b1500"],
  span_only: ["#c2b5d6", "#21173a"], strand_mismatch: ["#8c510a", "#fff"], missed: ["#b2182b", "#fff"] };
function renderOutcomes(runs) {
  const el = $("#outcomeChart");
  $("#outcomeContext").textContent = `(policy ${S.policy}, ${scopeInfo().label})`;
  const W = Math.min(1100, el.clientWidth || 1000), L = 230, R = 12, bh = 22, gap = 8;
  const x = (f) => L + f * (W - L - R);
  let y = 4, out = "";
  for (const r of runs) {
    const st = roleStyle(r), gate = resultGate(r), mid = y + bh / 2;
    out += marker(st.shape, 10, mid, st.shape === "ring" ? 4 : 5, r.color) + `<text x="22" y="${mid + 4}" font-weight="${st.rank < 2 ? 700 : 400}">${esc(runName(r).slice(0, 30))}</text>`;
    const rows = OUTCOME_ORDER.map((k) => [k, metric(r.run_id, "cds_outcome_" + k)]);
    if (!gate.ok) out += `<text x="${L}" y="${mid + 4}" class="lbl-muted" font-style="italic">${esc(gate.text)}${gate.why ? " — " + esc(gate.why.slice(0, 90)) : ""}</text>`;
    else if (rows.some(([, m]) => !m || m.status !== "ok")) out += `<text x="${L}" y="${mid + 4}" class="lbl-muted" font-style="italic">not available: the comparator reported no CDS outcome rows for this result${result(r.run_id)?.status === "legacy" ? " (pre-contract)" : ""}</text>`;
    else {
      let acc = 0;
      for (const [k, m] of rows) {
        const [fill, ink] = OUTCOME_FILL[k], w = x(acc + m.value) - x(acc);
        out += `<rect x="${x(acc).toFixed(1)}" y="${y}" width="${Math.max(0, w).toFixed(1)}" height="${bh}" fill="${fill}" data-tip="${esc(`${runName(r)} · ${OUTCOME[k][1]}\n${fmtInt(m.numerator)} of ${fmtInt(m.denominator)} reference genes (${fmtPct(m.value)})\nmetric cds_outcome_${k}`)}"/>`;
        if (w > 38) out += `<text x="${(x(acc) + w / 2).toFixed(1)}" y="${mid + 4}" text-anchor="middle" style="fill:${ink}" font-size="10.5" pointer-events="none">${fmtPct(m.value, 0)}</text>`;
        acc += m.value;
      }
    }
    y += bh + gap;
  }
  const legend = `<div class="legend">${OUTCOME_ORDER.map((k) => `<span><span class="swatch" style="background:${OUTCOME_FILL[k][0]}"></span> ${esc(OUTCOME[k][0])} <span class="small">(${esc(OUTCOME[k][1])})</span></span>`).join("")}</div>`;
  el.innerHTML = runs.length ? legend + svg(W, y, out, "reference genes by CDS outcome") : "";
}
function renderLimitBox(runs) {
  const byText = new Map();  // one entry per limitation, naming every run it applies to
  for (const r of runs) for (const l of runLimitations(r)) byText.set(l, [...(byText.get(l) || []), runName(r)]);
  const items = [...byText].map(([l, names]) => `<li><b>${esc(names.join(", "))}</b>: ${esc(l)}</li>`);
  for (const c of refMeta().caveats || []) items.push(`<li><b>Reference</b>: ${esc(c)}</li>`);
  $("#limitBox").innerHTML = items.length ? `<details class="panel limits" open><summary>Limitations relevant to this view (${items.length})</summary><ul>${items.join("")}</ul></details>` : "";
}

// ---------------------------------------------------------------- structural view
function renderStructure() {
  const runs = visibleRuns();
  const head = `<tr><th>Prediction run</th>${STRUCTURE.map((c) => `<th class="num"><span class="info" tabindex="0" data-tip="${esc(S.dict[c.id]?.definition || "")}">${esc(c.label)}</span><div class="cell-sub">${esc(c.sub)}</div></th>`).join("")}</tr>`;
  const body = runs.map((r) => `<tr class="${isSupp(r) ? "supp" : ""}">${runCell(r)}${STRUCTURE.map((c) => valueCell(r, c.id, "pct")).join("")}</tr>`).join("");
  $("#structureSummary").innerHTML = `<h2>Structural agreement <span class="small muted">(policy ${esc(S.policy)}, ${esc(scopeInfo().label)})</span></h2><div class="xscroll"><table class="grid summary">${head}${body}</table></div>`;
}
const LEN_SHORT = ["< 10", "10–30", "30–60", "60–100", "100–200", "200–500", "0.5–1k", "1–2k", "2–5k", "5–10k", "10–20k", "20–50k", "50–100k", "≥ 100k"];
// A segmented control: picking a button marks it selected (aria-checked, styled like the other active
// controls) before onPick runs, so the visible state always matches the state the plots use.
function seg(el, options, value, onPick) {
  el.innerHTML = options.map(([v, label, t]) => `<button type="button" role="radio" aria-checked="${v === value}" tabindex="${v === value ? 0 : -1}" data-v="${esc(v)}" ${t ? `data-tip="${esc(t)}"` : ""}>${esc(label)}</button>`).join("");
  el._pick = (v) => { $$("button", el).forEach((b) => { const on = b.dataset.v === v; b.setAttribute("aria-checked", on); b.tabIndex = on ? 0 : -1; }); onPick(v); };
  $$("button", el).forEach((b) => (b.onclick = () => el._pick(b.dataset.v)));
  if (!el.dataset.nav) { arrowNav(el, "button", (b) => el._pick(b.dataset.v)); el.dataset.nav = "1"; }
}
function spliceAvailability() {
  const rs = visibleRuns().map((r) => result(r.run_id)).filter(shown);
  const a = rs.map((r) => r.analyses?.splice_position_length).filter(Boolean);
  if (a.length && a.every((x) => x.status !== "ok")) return { ok: false, reason: a[0].reason };
  if (!a.length) return { ok: false, reason: rs.some((r) => r.status === "legacy") ? "pre_contract_result" : "no result" };
  return { ok: true };
}
async function renderPlots() {
  const p = policyInfo(), av = spliceAvailability();
  if (!av.ok) {
    const others = S.meta.policies.filter((x) => x.reference_transcript_selection !== "all");
    $("#plotStatus").innerHTML = `<div class="issue info">Not available under policy ${esc(p.id)} (${esc(av.reason)}). ${p.reference_transcript_selection === "all"
      ? "Splice-position and intron-length analyses need one reference transcript per gene; under policy A a gene with several isoforms would contribute shared introns several times, so the comparator does not compute them."
      : "The comparator did not report these analyses for these results."}
      ${others.length ? `<span>Show them for ${others.map((x) => `<button class="btn" type="button" data-pol="${esc(x.id)}">policy ${esc(x.id)}</button>`).join(" ")}</span>` : ""}</div>`;
    $$("#plotStatus [data-pol]").forEach((b) => (b.onclick = () => { S.policy = b.dataset.pol; renderAll(); }));
    $("#plotArea").hidden = true; return;
  }
  $("#plotArea").hidden = false; $("#plotStatus").innerHTML = "";
  seg($("#posCoord"), [["cds_span", "genomic CDS", "genomic CDS span: genomic bases from the coding start to the junction (introns included) over the span from coding start to coding end, in the direction of transcription"],
    ["coding", "coding sequence", "spliced coding sequence: CDS bases upstream of the junction over the CDS length (introns excluded); a donor and the acceptor of the same intron share one position"]], S.posCoord, (v) => { S.posCoord = v; drawPlots(); });
  seg($("#posMin"), [["0", "all CDS introns", "every reference CDS intron"], ["60", "introns ≥ 60 bp", "leave out reference CDS introns shorter than 60 bp"]], String(S.posMin), (v) => { S.posMin = Number(v); drawPlots(); });
  seg($("#lenMeasure"), [["exact", "complete intron", "both boundaries match"], ["donor", "donor site", "the 5' splice site matches"], ["acceptor", "acceptor site", "the 3' splice site matches"]], S.lenMeasure, (v) => { S.lenMeasure = v; drawPlots(); });
  const key = `${S.ref}|${S.policy}|${S.scope}`;
  if (!S.splice[key]) {
    try { S.splice[key] = await api("/api/splice", { reference: S.ref, policy: S.policy, scope: S.scope }); }
    catch (e) { S.splice[key] = { error: e.message }; }
  }
  drawPlots();
}
function drawPlots() {
  const data = S.splice[`${S.ref}|${S.policy}|${S.scope}`];
  if (!data) return;
  for (const id of ["#positionChart", "#positionTable", "#lengthChart", "#lengthTable", "#positionNote", "#lengthNote", "#spliceDownloads"]) $(id).innerHTML = "";
  if (data.error) { $("#plotStatus").innerHTML = `<p class="na">Not available: ${esc(data.error)}</p>`; return; }
  const runs = visibleRuns().filter((r) => data.sites.some((x) => x.run_id === r.run_id));
  const missing = visibleRuns().filter((r) => !runs.includes(r)).map((r) => `${runName(r)}: ${result(r.run_id)?.analyses?.splice_position_length?.reason || resultGate(r).text}`);
  $("#plotStatus").innerHTML = missing.length ? `<div class="issue warning">No splice tables for: ${missing.map(esc).join("; ")}</div>` : "";
  if (!runs.length) return;
  drawPositionChart(data, runs); drawLengthChart(data, runs); writeHash();
  $("#spliceDownloads").innerHTML = `<span class="small muted">Comparator tables (TSV, policy ${esc(S.policy)}, ${esc(scopeInfo().label)}):</span>
    <a class="btn" href="${esc(apiUrl("/api/splice", { reference: S.ref, policy: S.policy, scope: S.scope, format: "csv", table: "sites" }))}">splice sites by position</a>
    <a class="btn" href="${esc(apiUrl("/api/splice", { reference: S.ref, policy: S.policy, scope: S.scope, format: "csv", table: "lengths" }))}">introns by length</a>`;
}
function niceMax(v) { return Math.min(1, Math.max(0.1, Math.ceil(v * 10 + 1e-9) / 10)); }
function svg(w, h, inner, label = "") { return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img" aria-label="${esc(label)}">${inner}</svg>`; }
function lineChart({ width, height, labels, series, yMax, yTitle, xTitle, footer, edgeLabels, label }) {
  const L = 52, R = 14, T = 14, B = 40 + (footer ? 18 : 0), n = labels.length;
  const pw = width - L - R, ph = height - T - B;
  const x = (i) => L + ((i + 0.5) * pw) / n, y = (v) => T + ph - (v / yMax) * ph;
  let out = `<g class="grid">`;
  const step = yMax > 0.5 ? 0.2 : yMax > 0.2 ? 0.1 : 0.05;
  for (let t = 0; t <= yMax + 1e-9; t += step) out += `<line x1="${L}" x2="${width - R}" y1="${y(t)}" y2="${y(t)}"/><text x="${L - 6}" y="${y(t) + 4}" text-anchor="end" class="lbl-muted">${Math.round(t * 100)}%</text>`;
  out += `</g><line x1="${L}" x2="${width - R}" y1="${T + ph}" y2="${T + ph}" stroke="var(--line-strong)"/>`;
  if (edgeLabels) for (let i = 0; i <= n; i += 2) out += `<text x="${L + (i * pw) / n}" y="${T + ph + 16}" text-anchor="middle" class="lbl-muted">${edgeLabels(i)}</text>`;
  else labels.forEach((l, i) => (out += `<text x="${x(i)}" y="${T + ph + 16}" text-anchor="middle" class="lbl-muted">${esc(l)}</text>`));
  if (footer) footer.forEach((f, i) => (out += `<text x="${x(i)}" y="${T + ph + 32}" text-anchor="middle" class="lbl-muted" font-size="10">${esc(f)}</text>`));
  out += `<text x="${L + pw / 2}" y="${height - 4}" text-anchor="middle" class="lbl-muted">${esc(xTitle)}</text>`;
  out += `<text transform="translate(12 ${T + ph / 2}) rotate(-90)" text-anchor="middle" class="lbl-muted">${esc(yTitle)}</text>`;
  // Others first, then the primary comparator, then the focus run, so the focus is never hidden.
  for (const s of [...series].sort((a, b) => (b.style?.rank ?? 2) - (a.style?.rank ?? 2))) {
    const st = s.style || OTHER_STYLE;
    let path = "", pen = false;
    s.points.forEach((p, i) => { if (p.v === null) { pen = false; return; } path += `${pen ? "L" : "M"}${x(i).toFixed(1)},${y(p.v).toFixed(1)}`; pen = true; });
    if (st.rank < 2) out += `<path d="${path}" fill="none" stroke="#fff" stroke-width="${st.width + 3}" stroke-linejoin="round" opacity="0.9"/>`;
    out += `<path d="${path}" fill="none" stroke="${s.color}" stroke-width="${st.width}" stroke-linejoin="round" ${s.dash ? 'stroke-dasharray="5 3"' : ""}/>`;
    s.points.forEach((p, i) => { if (p.v !== null) out += marker(st.shape, x(i), y(p.v), st.r, s.color, p.tip); });
  }
  return svg(width, height, out, label);
}
function legendHtml(runs) { return `<div class="legend">${runs.map((r) => `<span class="${ROLE_STYLE[r.display_role] ? "legend-role" : ""}">${swatch(r)} ${esc(runName(r))}${roleTag(r) ? " " + roleTag(r) : ""}${isSupp(r) ? " (supplementary, dashed)" : ""}</span>`).join("")}</div>`; }
function rate(num, den) { return den ? num / den : null; }
function drawPositionChart(data, runs) {
  const el = $("#positionChart");
  const W = Math.min(1100, el.clientWidth || 1000), half = W >= 760 ? Math.floor((W - 16) / 2) : W;
  const labels = Array.from({ length: 10 }, (_, i) => `${(i / 10).toFixed(1)}–${((i + 1) / 10).toFixed(1)}`);
  const rows = (run, kind) => data.sites.filter((x) => x.run_id === run && (x.min_intron_length ?? 0) === S.posMin && x.coordinate === S.posCoord && x.kind === kind).sort((a, b) => a.bin - b.bin);
  let max = 0;
  for (const r of runs) for (const k of ["donor", "acceptor"]) for (const b of rows(r.run_id, k)) { const v = rate(b.assessed - b.supported, b.assessed); if (v !== null) max = Math.max(max, v); }
  const yMax = niceMax(max);
  const panel = (kind) => lineChart({ width: half, height: 250, labels, yMax, label: `unsupported ${kind} sites by position`, yTitle: "unsupported (% of assessed)",
    xTitle: `${kind} sites · ${S.posCoord === "coding" ? "position in the spliced coding sequence (CDS bases only)" : "position on the genomic CDS span (introns included)"}`, edgeLabels: (i) => (i / 10).toFixed(1),
    series: runs.map((r) => ({ color: r.color, dash: isSupp(r), style: roleStyle(r), points: rows(r.run_id, kind).map((b) => ({ v: rate(b.assessed - b.supported, b.assessed),
      tip: `${runName(r)} · ${kind} sites, position ${labels[b.bin]}\nunsupported ${fmtInt(b.assessed - b.supported)} of ${fmtInt(b.assessed)} assessed\n${fmtInt(b.missing_locus)} further sites in genes whose locus was not recovered (excluded)` })) })) });
  el.innerHTML = legendHtml(runs) + `<div class="panels">${panel("donor")}${panel("acceptor")}</div>`;
  $("#positionNote").innerHTML = `${S.posCoord === "coding" ? "<b>Coding sequence</b>: CDS bases upstream of the junction ÷ CDS length, in the direction of transcription (introns excluded; a donor and its acceptor share one position)." : "<b>Genomic CDS</b>: genomic bases from the coding start to the junction ÷ genomic span from coding start to coding end, in the direction of transcription (introns included, so long introns move later sites towards 1)."} 0 = coding start, 1 = coding end; both are computed by the comparator for one reference transcript per gene. A reference donor (acceptor) is <b>unsupported</b> when no evaluated predicted transcript on the same sequence and strand has a CDS intron with that coordinate (comparator support rule). Rates exclude sites of genes whose locus was not recovered (missing-locus sites, counted in the table). Empty bins are not drawn. ${S.posMin ? `Introns shorter than ${S.posMin} bp are left out.` : ""}`;
  const head = `<tr><th>Prediction run</th><th>site</th>${labels.map((l) => `<th class="num">${l}</th>`).join("")}<th class="num">missing-locus sites</th></tr>`;
  const body = runs.flatMap((r) => ["donor", "acceptor"].map((k) => { const b = rows(r.run_id, k);
    return `<tr><td>${esc(runName(r))}</td><td>${k}</td>${b.map((x) => x.assessed ? `<td class="num">${fmtPct((x.assessed - x.supported) / x.assessed)}<div class="cell-sub">${fmtInt(x.assessed - x.supported)} / ${fmtInt(x.assessed)}</div></td>` : '<td class="num"><span class="na">not available</span><div class="cell-sub">0 assessed</div></td>').join("")}
      <td class="num">${fmtInt(b.reduce((t, x) => t + x.missing_locus, 0))}</td></tr>`; })).join("");
  $("#positionTable").innerHTML = `<table class="grid compact">${head}${body}</table>`;
}
function drawLengthChart(data, runs) {
  const el = $("#lengthChart");
  const W = Math.min(1100, el.clientWidth || 1000);
  const rows = (run) => data.lengths.filter((x) => x.run_id === run).sort((a, b) => a.bin - b.bin);
  const first = rows(runs[0].run_id);
  const refTotals = first.map((b) => b.eligible + b.missing_locus);
  const m = S.lenMeasure, word = { exact: "complete introns recovered", donor: "donor sites supported", acceptor: "acceptor sites supported" }[m];
  el.innerHTML = legendHtml(runs) + lineChart({ width: W, height: 290, labels: first.map((b) => LEN_SHORT[b.bin] || b.label), yMax: 1, label: word, yTitle: `${word} (% of eligible)`,
    xTitle: "reference CDS intron length (bp, fixed bins)", footer: refTotals.map((n) => `n=${fmtInt(n)}`),
    series: runs.map((r) => ({ color: r.color, dash: isSupp(r), style: roleStyle(r), points: rows(r.run_id).map((b) => ({ v: rate(b[m], b.eligible),
      tip: `${runName(r)} · introns ${b.label}\n${word}: ${fmtInt(b[m])} of ${fmtInt(b.eligible)}\n${fmtInt(b.missing_locus)} introns in genes whose locus was not recovered (excluded)` })) })) });
  $("#lengthNote").innerHTML = "n = reference CDS introns per bin (one transcript per gene). Associations, not causes: bins hold different genes and different numbers of introns.";
  const head = `<tr><th>Prediction run</th>${first.map((b) => `<th class="num">${esc(b.label)}</th>`).join("")}</tr>`;
  const body = runs.map((r) => `<tr><td>${esc(runName(r))}</td>${rows(r.run_id).map((b) => b.eligible ? `<td class="num">${fmtPct(b[m] / b.eligible)}<div class="cell-sub">${fmtInt(b[m])} / ${fmtInt(b.eligible)}</div></td>` : `<td class="num"><span class="na">not available</span><div class="cell-sub">0 eligible</div></td>`).join("")}</tr>`).join("");
  $("#lengthTable").innerHTML = `<table class="grid compact">${head}${body}</table>`;
}

// ---------------------------------------------------------------- all metrics
function renderDetail() {
  const runs = visibleRuns();
  const groups = S.meta.metric_dictionary.groups;
  const byGroup = {};
  for (const m of S.meta.metric_dictionary.metrics) (byGroup[m.group] = byGroup[m.group] || []).push(m);
  const defs = S.meta.metric_definitions ? Object.values(S.meta.metric_definitions)[0] || {} : {};
  $("#detailTables").innerHTML = groups.map((g) => {
    const ms = byGroup[g.id] || [];
    if (!ms.length) return "";
    const head = `<tr><th>Metric</th><th>Denominator</th>${runs.map((r) => `<th class="num">${esc(runName(r))}${isSupp(r) ? '<div class="cell-sub">supplementary</div>' : ""}</th>`).join("")}</tr>`;
    const body = ms.map((d) => {
      const def = defs[d.id] || {};
      const t = [d.definition, def.numerator ? "Numerator: " + def.numerator : "", def.denominator_id ? "Denominator: " + def.denominator_id : "", ...(d.limitations || []).map((l) => "Limitation: " + l)].filter(Boolean).join("\n");
      return `<tr><td><span class="info" tabindex="0" data-tip="${esc(t)}">${esc(d.label)}</span><div class="cell-sub mono">${esc(d.id)}</div></td><td class="small muted">${esc(def.denominator_id || (d.group === "historical" ? "R / Q" : "—"))}</td>
        ${runs.map((r) => valueCell(r, d.id, /_f1(_|$)/.test(d.id) ? "f1" : (metric(r.run_id, d.id)?.denominator ?? null) === null ? "count" : "pct", { bar: false })).join("")}</tr>`;
    }).join("");
    return `<h3>${esc(g.label)}</h3><div class="tablewrap"><table class="grid compact">${head}${body}</table></div>`;
  }).join("");
  $("#dlMetricsView").href = apiUrl("/api/overview", { reference: S.ref, format: "csv", policy: S.policy, scope: S.scope });
  $("#dlMetricsAll").href = apiUrl("/api/overview", { reference: S.ref, format: "csv" });
  const sel = $("#policyMetric");
  if (!sel.options.length) {
    sel.innerHTML = [...HEADLINE, ...STRUCTURE].map((c) => `<option value="${c.id}">${esc(S.dict[c.id]?.label || c.label)}</option>`).join("");
    sel.value = "exact_cds_one_to_one_f1"; sel.onchange = renderPolicyChart;
  }
  renderPolicyChart(); renderDownloads();
}
function renderPolicyChart() {
  const id = $("#policyMetric").value, runs = visibleRuns(), el = $("#policyChart");
  const W = Math.min(980, el.clientWidth || 980), L = 260, bh = 16;
  let y = 8, out = "";
  const x = (v) => L + v * (W - L - 150);
  for (const p of S.meta.policies) {
    out += `<text x="4" y="${y + 12}" font-weight="700">${esc(p.id)} · ${esc(p.label)}</text>`; y += 18;
    for (const r of runs) {
      const res = result(r.run_id, p.id), m = metric(r.run_id, id, p.id);
      const st = roleStyle(r);
      out += marker(st.shape, 20, y + 8, st.shape === "ring" ? 3.5 : 4.5, r.color) + `<text x="30" y="${y + 12}" ${st.rank < 2 ? 'font-weight="700"' : 'class="lbl-muted"'}>${esc(runName(r).slice(0, 34))}</text>`;
      if (!shown(res)) out += `<text x="${L}" y="${y + 12}" class="lbl-muted" font-style="italic">${esc(res ? RESULT_LABEL[res.status] : "no result")}</text>`;
      else if (!m || m.status !== "ok") out += `<text x="${L}" y="${y + 12}" class="lbl-muted" font-style="italic">${esc(m ? `${NA_TEXT[m.status] || m.status} (${m.reason})` : "not available")}</text>`;
      else {
        out += `<rect x="${L}" y="${y + 2}" width="${Math.max(1, x(m.value) - L)}" height="${bh - 4}" fill="${r.color}" data-tip="${esc(metricTip(m, id))}"></rect>`;
        out += `<text x="${x(m.value) + 6}" y="${y + 12}">${id.endsWith("_f1") ? m.value.toFixed(3) : fmtPct(m.value)} <tspan class="lbl-muted">(${fmtInt(m.numerator)} / ${fmtInt(m.denominator)})</tspan></text>`;
      }
      y += bh + 4;
    }
    y += 12;
  }
  el.innerHTML = svg(W, y, out, "policy effect");
}
function renderDownloads() {
  const runs = visibleRuns().filter((r) => shown(result(r.run_id)));
  $("#downloads").innerHTML = runs.length ? `<table class="grid compact"><tr><th>Prediction run</th><th>Downloads (policy ${esc(S.policy)})</th></tr>${runs.map((r) => {
    const res = result(r.run_id), chk = res.pairs_check || {};
    const links = [`<a class="btn" href="${esc(apiUrl("/api/pairs", { reference: S.ref, policy: S.policy, run: r.run_id, format: "csv" }))}" data-tip="${esc(`${chk.rows ?? "?"} rows; TP ${chk.true_positives ?? "?"}; check ${chk.status || "?"}`)}">1:1 exact pairs</a>`,
      `<a class="btn" href="${esc(apiUrl("/api/per_sequence", { reference: S.ref, policy: S.policy, run: r.run_id, format: "csv" }))}">per sequence</a>`];
    if (harmonised(r)) links.push(`<a class="btn" href="${esc(apiUrl("/api/stop_audit", { run: r.run_id, format: "csv" }))}">stop-codon audit</a>`);
    return `<tr>${runCell(r)}<td>${links.join(" ")}</td></tr>`;
  }).join("")}</table><p class="muted small">Per-sequence counts are drill-down only; scope values come from the comparator's own scope rows, never from summing sequences.</p>` : '<p class="na">No displayed results.</p>';
}

// ---------------------------------------------------------------- genes
function setupGeneFilters() {
  const runs = visibleRuns();
  $("#focusRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(runName(r))}</option>`).join("");
  $("#compareRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(runName(r))}</option>`).join("");
  if (runs[1]) $("#compareRun").value = runs[1].run_id;
  $("#statusChecks").innerHTML = '<legend class="sr">CDS outcome in the selected run</legend>' + Object.entries(OUTCOME).map(([k, [lab, t]]) => `<label data-tip="${esc(t)}"><input type="checkbox" value="${k}"> <span class="status st-${k}">${esc(lab)}</span></label>`).join("");
  api("/api/chromosomes", { reference: S.ref }).then((chroms) => {
    $("#chrom").innerHTML = '<option value="">all sequences</option>' + chroms.map((c) => `<option>${esc(c)}</option>`).join("");
  });
  $("#geneFilters").onsubmit = (e) => { e.preventDefault(); S.genePage = 0; loadGenes(); };
  $("#resetFilters").onclick = () => { $("#geneFilters").reset(); S.genePage = 0; loadGenes(); };
  $("#dlGenes").onclick = () => (location.href = apiUrl("/api/genes", { ...geneParams(), format: "csv", runs: visibleRuns().map((r) => r.run_id).join(",") }));
}
function geneParams() {
  return { reference: S.ref, policy: S.policy, scope: S.scope, q: $("#q").value, chrom: $("#chrom").value, multi: $("#multi").value,
    focus_run: $("#focusRun").value, status: $$("#statusChecks input:checked").map((c) => c.value).join(","), split: $("#splitOnly").checked ? "1" : "",
    one_to_one: $("#oneToOne").value, contrast: $("#contrast").value, compare_run: $("#compareRun").value, limit: S.pageSize, offset: S.genePage * S.pageSize };
}
async function loadGenes() {
  S.genesLoaded = true;
  const token = (S.geneToken = (S.geneToken || 0) + 1);
  try { const data = await api("/api/genes", geneParams()); if (token !== S.geneToken) return; S.lastGenes = data; renderGeneTable(data); }
  catch (e) { showError($("#geneTable"), e); }
}
function outcomeChip(o, run) {
  const res = result(run.run_id);
  if (!o) return `<span class="status st-none" data-tip="${esc(shown(res) ? "gene not in this result" : (RESULT_LABEL[res?.status] || "no result"))}">—</span>`;
  if (!o.cds_outcome) return `<span class="status st-none" data-tip="${esc(`pre-contract result: no comparator CDS outcome category\ncomparator class ${o.classification} / CDS ${o.classification_cds}\ncoordinate-exact (best pair): ${o.cds_exact ? "yes" : "no"}`)}">${o.cds_exact ? "exact*" : o.classification}</span>`;
  const [lab, t] = OUTCOME[o.cds_outcome] || [o.cds_outcome, ""];
  const extra = [`comparator: ${o.classification} / CDS ${o.classification_cds}`, `shared CDS bases ${fmtInt(o.cds_shared_bases)} (reciprocal overlap ${o.cds_overlap})`,
    o.cds_matched_id ? `best CDS partner ${o.cds_matched_id}` : "", o.one_to_one_partner ? `one-to-one exact pair: ${o.one_to_one_partner}` : o.cds_exact ? "coordinate-exact, but not in the one-to-one matching" : "",
    o.counterpart_count >= 2 ? `split: ${o.counterpart_ids}` : "", `CDS intron chain: best pair ${o.cds_intron_chain_match}, any pair ${o.cds_intron_chain_any_pair ?? "NA"}`].filter(Boolean).join("\n");
  return `<span class="status st-${esc(o.cds_outcome)}" data-tip="${esc(t + "\n" + extra)}">${esc(lab)}${o.counterpart_count >= 2 ? " ⇉" : ""}${o.cds_exact && !o.one_to_one_partner ? " ≠1:1" : ""}</span>`;
}
function renderGeneTable(data) {
  const runs = visibleRuns();
  $("#geneCount").textContent = `${fmtInt(data.total)} evaluated reference genes match (policy ${S.policy}, ${scopeInfo().label}); showing ${data.total ? data.offset + 1 : 0}–${data.offset + data.rows.length}. ⇉ split · ≠1:1 exact but outside the one-to-one matching. Enter or click opens the locus.`;
  const head = `<tr><th>Reference gene</th><th>Name</th><th>Location</th><th>CDS</th>${runs.map((r) => `<th>${swatch(r)} ${esc(runName(r))}</th>`).join("")}</tr>`;
  const body = data.rows.map((g) => `<tr class="click ${S.selectedGene === g.gene_id ? "selected" : ""}" tabindex="0" data-gene="${esc(g.gene_id)}">
      <td class="mono">${esc(g.gene_id)}${g.original_gene_id && g.original_gene_id !== g.gene_id ? `<div class="cell-sub">file ID ${esc(g.original_gene_id)}</div>` : ""}</td>
      <td>${esc(g.name || "")}</td><td class="mono">${esc(g.chrom)}:${fmtInt(g.start)}-${fmtInt(g.end)} ${esc(g.strand)}</td>
      <td class="small">${g.max_cds_count === null ? "—" : (g.max_cds_count > 1 ? "multi" : "single") + "-segment"} · ${g.evaluated_transcript_count ?? "?"} evaluated</td>
      ${runs.map((r) => `<td>${outcomeChip(g.outcomes[r.run_id], r)}</td>`).join("")}</tr>`).join("");
  $("#geneTable").innerHTML = `<table class="grid">${head}${body}</table>`;
  $$("#geneTable tr.click").forEach((tr) => {
    const open = () => { const g = data.rows.find((x) => x.gene_id === tr.dataset.gene); S.selectedGene = g.gene_id;
      $$("#geneTable tr").forEach((x) => x.classList.toggle("selected", x === tr)); openLocus(g.chrom, g.start, g.end, g.gene_id); };
    tr.onclick = open; tr.onkeydown = (e) => { if (e.key === "Enter") open(); };
  });
  pager($("#genePager"), data.total, S.genePage, (p) => { S.genePage = p; loadGenes(); });
}
function pager(el, total, page, go) {
  const pages = Math.max(1, Math.ceil(total / S.pageSize));
  el.innerHTML = `<button class="btn" ${page <= 0 ? "disabled" : ""} data-p="${page - 1}">‹ previous</button><span class="small">page ${page + 1} of ${pages}</span>
    <button class="btn" ${page >= pages - 1 ? "disabled" : ""} data-p="${page + 1}">next ›</button>`;
  $$("button", el).forEach((b) => (b.onclick = () => go(Number(b.dataset.p))));
}

// ---------------------------------------------------------------- predictions
function setupPredFilters() {
  $("#predRun").innerHTML = visibleRuns().map((r) => `<option value="${esc(r.run_id)}">${esc(runName(r))}</option>`).join("");
  $("#pnovel").innerHTML = '<option value="">any Novel context</option>' + NOVEL.map((id) => `<option value="${id.replace("query_novel_", "")}">${esc(S.dict[id]?.label || id)}</option>`).join("");
  $("#predFilters").onsubmit = (e) => { e.preventDefault(); S.predPage = 0; loadPredictions(); };
  $("#dlPreds").onclick = () => (location.href = apiUrl("/api/predictions", { ...predParams(), format: "csv" }));
}
function predParams() {
  return { reference: S.ref, policy: S.policy, scope: S.scope, run: $("#predRun").value, q: $("#pq").value, classification: $("#pclass").value,
    novel_category: $("#pnovel").value, exact: $("#pexact").value, one_to_one: $("#pone").value, merge: $("#pmerge").checked ? "1" : "", limit: S.pageSize, offset: S.predPage * S.pageSize };
}
async function loadPredictions() {
  S.predsLoaded = true;
  const run = runById($("#predRun").value), res = run && result(run.run_id);
  if (!shown(res)) { $("#predTable").innerHTML = `<p class="na">No displayed result for this run and policy (${esc(res ? RESULT_LABEL[res.status] : "no result")}).</p>`; $("#predCount").textContent = ""; $("#predPager").innerHTML = ""; return; }
  const token = (S.predToken = (S.predToken || 0) + 1);
  try {
    const data = await api("/api/predictions", predParams());
    if (token !== S.predToken) return;
    $("#predCount").textContent = `${fmtInt(data.total)} predicted genes match; showing ${data.total ? data.offset + 1 : 0}–${data.offset + data.rows.length}. Novel means unmatched relative to the evaluated reference, not a demonstrated false positive.`;
    const head = "<tr><th>Predicted gene (comparison ID)</th><th>Location</th><th>Classification</th><th>Best-pair exact</th><th>1:1 exact pair</th><th>Reference partner</th><th>Novel context</th><th>Counterparts</th></tr>";
    const body = data.rows.map((q) => `<tr class="click" tabindex="0" data-chrom="${esc(q.chrom)}" data-s="${q.start}" data-e="${q.end}">
      <td class="mono">${esc(q.gene_id)}${q.original_gene_id && q.original_gene_id !== q.gene_id ? `<div class="cell-sub">file ID ${esc(q.original_gene_id)}</div>` : ""}</td>
      <td class="mono">${esc(q.chrom)}:${fmtInt(q.start)}-${fmtInt(q.end)} ${esc(q.strand)}</td>
      <td>${esc(q.classification)}<div class="cell-sub">CDS ${esc(q.classification_cds)}</div></td>
      <td>${q.cds_exact ? '<span class="status st-exact">exact</span>' : '<span class="small muted">no</span>'}</td>
      <td class="mono small">${esc(q.one_to_one_partner || "—")}</td>
      <td class="mono">${esc(q.cds_matched_id || "")}<div class="cell-sub">${esc(q.best_cds_ref_tx || "")}</div></td>
      <td class="small">${esc((q.novel_category || "").replaceAll("_", " "))}</td>
      <td class="small">${q.counterpart_count >= 2 ? `<span class="tag warn">merge</span> ${esc(q.counterpart_ids)}` : esc(q.counterpart_ids || "")}</td></tr>`).join("");
    $("#predTable").innerHTML = `<table class="grid">${head}${body}</table>`;
    $$("#predTable tr.click").forEach((tr) => { const open = () => { selectTab("genes"); openLocus(tr.dataset.chrom, Number(tr.dataset.s), Number(tr.dataset.e), null); };
      tr.onclick = open; tr.onkeydown = (e) => { if (e.key === "Enter") open(); }; });
    pager($("#predPager"), data.total, S.predPage, (p) => { S.predPage = p; loadPredictions(); });
  } catch (e) { showError($("#predTable"), e); }
}

// ---------------------------------------------------------------- locus
function openLocus(chrom, start, end, geneId) {
  $("#locusPanel").hidden = false;
  const pad = Math.max(500, Math.round((end - start) * 0.15));
  S.locus = { chrom, start: Math.max(1, start - pad), end: end + pad, gene: geneId };
  drawLocus();
  $("#locusPanel").scrollIntoView({ behavior: "smooth", block: "start" });
  $("#locusTitle").focus({ preventScroll: true });
}
$("#locusInput").addEventListener("keydown", (e) => {
  if (e.key !== "Enter") return;
  const m = /^\s*([^:\s]+):([\d,]+)-([\d,]+)\s*$/.exec(e.target.value);
  if (!m) { $("#locusNotes").textContent = "Use seq:start-end"; return; }
  S.locus = { chrom: m[1], start: Number(m[2].replaceAll(",", "")), end: Number(m[3].replaceAll(",", "")), gene: S.locus?.gene };
  drawLocus();
});
$$("[data-zoom]").forEach((b) => (b.onclick = () => { if (!S.locus) return; const f = Number(b.dataset.zoom); const mid = (S.locus.start + S.locus.end) / 2;
  const half = Math.max(50, ((S.locus.end - S.locus.start) * f) / 2); S.locus.start = Math.max(1, Math.round(mid - half)); S.locus.end = Math.round(mid + half); drawLocus(); }));
$$("[data-move]").forEach((b) => (b.onclick = () => { if (!S.locus) return; const d = Math.round((S.locus.end - S.locus.start) * Number(b.dataset.move));
  S.locus.start = Math.max(1, S.locus.start + d); S.locus.end += d; drawLocus(); }));
$("#isoforms").onchange = () => S.locus && drawLocus();
function ivs(text) { return text ? text.split(",").map((p) => p.split("-").map(Number)) : []; }
function pack(items, key) {
  const rows = [];
  for (const it of items) { const [s, e] = key(it); let row = rows.findIndex((last) => last < s - 4); if (row < 0) { rows.push(e); row = rows.length - 1; } else rows[row] = e; it._row = row; }
  return rows.length;
}
function excludedText(t) { return t.evaluated ? `evaluated under policy ${S.policy} (selected by ${t.selected_by || "?"})` : `not evaluated under policy ${S.policy}: ${t.exclusion_reason || "excluded"}`; }
async function drawLocus() {
  const L = S.locus, runs = visibleRuns(), ref = refMeta();
  $("#locusInput").value = `${L.chrom}:${L.start}-${L.end}`;
  let data;
  try { data = await api("/api/locus", { reference: S.ref, chrom: L.chrom, start: L.start, end: L.end, policy: S.policy, runs: runs.map((r) => r.run_id).join(","), gene_id: L.gene || "", isoforms: $("#isoforms").value }); }
  catch (e) { showError($("#locusSvg"), e); return; }
  L.start = data.start; L.end = data.end;
  $("#locusTitle").textContent = `Locus ${data.chrom}:${fmtInt(data.start)}-${fmtInt(data.end)} (${fmtInt(data.end - data.start + 1)} bp) · policy ${S.policy}`;
  const notes = [];
  if (data.window_clamped) notes.push(`window limited to ${fmtInt(data.limits.max_window)} bp`);
  if (data.reference_transcripts_truncated) notes.push(`reference transcripts limited to ${data.limits.track_limit}`);
  if (data.reference_transcripts_hidden) notes.push(`${data.reference_transcripts_hidden} reference isoforms hidden (isoform selector)`);
  for (const [run, t] of Object.entries(data.runs)) { if (t.truncated) notes.push(`${runName(runById(run))}: limited to ${data.limits.track_limit} transcripts`); if (t.hidden_excluded) notes.push(`${runName(runById(run))}: ${t.hidden_excluded} excluded isoforms hidden`); }
  $("#locusNotes").textContent = notes.join(" · ");
  $("#locusLegend").innerHTML = `<span>Exact evaluated models rebuilt from the compared input files and the comparator's evaluated-transcript export (no re-selection).</span>
    <span>UTR exon thin outline · CDS thick · ★ Ensembl_canonical · ◆ compared isoform (best CDS pair) · faded, dashed = excluded under this policy · red tick = stop codon added by --add-stop-codon (audited)</span>`;
  const W = Math.max(760, $("#locusSvg").clientWidth || 1100), LBL = 250, R = W - 12, rowH = 20;
  const x = (p) => LBL + ((p - data.start) / (data.end - data.start + 1)) * (R - LBL);
  const clip = (p) => Math.min(R, Math.max(LBL, x(p)));
  let y = 6, out = "";
  const span = data.end - data.start + 1, step = Math.pow(10, Math.floor(Math.log10(span / 5)));
  const tickStep = span / step > 10 ? step * 2 : step;
  out += `<g class="axis"><line x1="${LBL}" y1="${y + 16}" x2="${R}" y2="${y + 16}"/>`;
  for (let t = Math.ceil(data.start / tickStep) * tickStep; t <= data.end; t += tickStep) out += `<line x1="${x(t)}" y1="${y + 12}" x2="${x(t)}" y2="${y + 20}"/><text x="${x(t)}" y="${y + 9}" text-anchor="middle" class="lbl-muted">${fmtInt(t)}</text>`;
  out += `</g>`; y += 30;
  const bandTop = y, focusGene = data.reference_genes.find((g) => g.gene_id === L.gene);
  if (!ref.models_available) {
    out += `<text x="4" y="${y + 12}" font-weight="700">Reference</text><text x="${LBL}" y="${y + 12}" class="lbl-muted">models unavailable: ${esc((ref.models_error || "").slice(0, 120))}</text>`; y += rowH;
  } else {
    const ctx = data.reference_genes;
    const nrow = pack(ctx, (g) => [x(g.start), x(g.start) + Math.max(x(g.end) - x(g.start), 7 * (g.name || g.gene_id).length)]);
    out += `<text x="4" y="${y + 12}" font-weight="700">Reference genes</text>`;
    for (const g of ctx) {
      const gy = y + g._row * 16;
      out += `<g data-tip="${esc(`${g.gene_id}${g.name ? " (" + g.name + ")" : ""}\n${g.biotype} ${g.chrom}:${g.start}-${g.end} ${g.strand}\n${g.evaluated ? "evaluated under at least one policy" : "not evaluated (biotype or evaluation mode)"}`)}">
        <rect x="${clip(g.start)}" y="${gy + 4}" width="${Math.max(2, clip(g.end) - clip(g.start))}" height="6" fill="${g.evaluated ? "var(--ref)" : "var(--line-strong)"}"/>
        <text x="${clip(g.start)}" y="${gy + 3}" font-size="10" class="${g.evaluated ? "" : "lbl-muted"}">${esc(g.name || g.gene_id)} ${g.strand === "-" ? "◀" : "▶"}</text></g>`;
    }
    y += Math.max(1, nrow) * 16 + 12;
    out += `<text x="4" y="${y + 12}" font-weight="700">Reference transcripts</text>`; y += 18;
    for (const t of data.reference_transcripts) {
      out += glyph(t, y, clip, "var(--ref)", !t.evaluated, `${t.canonical ? "★ " : ""}${t.partner_of_runs.length ? "◆ " : ""}${t.transcript_id}`,
        `${t.transcript_id} (${t.biotype || "?"})\ngene ${t.gene_id}\n${excludedText(t)}${t.tags ? "\ntags: " + t.tags : ""}${t.partner_of_runs.length ? "\nbest CDS pair for: " + t.partner_of_runs.map((r) => runName(runById(r))).join(", ") : ""}`,
        t.partner_of_runs.map((r) => runById(r)?.color), null);
      y += rowH;
    }
    if (!data.reference_transcripts.length) { out += `<text x="${LBL}" y="${y + 12}" class="lbl-muted">no reference transcripts in this window for the isoform selection</text>`; y += rowH; }
  }
  for (const r of runs) {
    const tr = data.runs[r.run_id];
    y += 6;
    out += `<rect x="4" y="${y + 3}" width="10" height="10" fill="${r.color}"/><text x="18" y="${y + 12}" font-weight="700">${esc(runName(r))} <tspan class="lbl-muted" font-weight="400">${harmonised(r) ? "terminal stop codons added by the comparator in red (audited)" : ""}</tspan></text>`;
    y += 18;
    const err = S.meta.model_errors?.["run:" + r.run_id];
    if (err || !tr || !tr.transcripts.length) {
      out += `<text x="${LBL}" y="${y + 12}" class="lbl-muted">${esc(err ? "models unavailable: " + err.slice(0, 110) : tr && tr.hidden_excluded ? "only excluded isoforms here (isoform selector)" : "no predictions in this window")}</text>`; y += rowH; continue;
    }
    for (const t of tr.transcripts) {
      const o = t.outcome, st = t.stop_codon;
      const lab = `${t.gene_id}${o?.one_to_one_partner ? " · 1:1 " + o.one_to_one_partner : o?.cds_exact ? " · exact" : ""}${o ? " · " + o.classification : ""}`;
      const stopTip = st ? `\nstop codon (--add-stop-codon): ${st.status} — ${st.reason}\nterminal CDS ${st.terminal_cds_start}-${st.terminal_cds_end}${st.status === "extended" ? " → " + st.new_terminal_cds_start + "-" + st.new_terminal_cds_end : ""}` : "";
      out += glyph(t, y, clip, r.color, !t.evaluated, lab, `${t.transcript_id}${t.original_transcript_id && t.original_transcript_id !== t.transcript_id ? " (file ID " + t.original_transcript_id + ")" : ""}\ngene ${t.gene_id}\n${excludedText(t)}\n${o ? `comparator: ${o.classification} / CDS ${o.classification_cds}${o.novel_category ? "\nNovel: " + o.novel_category : ""}${o.counterpart_count >= 2 ? "\nmerge (" + o.counterpart_count + " reference counterparts)" : ""}` : "no outcome"}${stopTip}`, [], st && st.status === "extended" ? st : null);
      y += rowH;
    }
  }
  if (focusGene) out = `<rect x="${clip(focusGene.start)}" y="${bandTop}" width="${Math.max(2, clip(focusGene.end) - clip(focusGene.start))}" height="${y - bandTop}" fill="var(--select)" opacity="0.7"/>` + out;
  $("#locusSvg").innerHTML = svg(W, y + 8, out, `locus ${data.chrom}:${data.start}-${data.end}`);
  if (L.gene) loadGeneDetail(L.gene); else $("#geneDetail").innerHTML = "";
  function glyph(t, gy, cl, color, excluded, label, tipText, markers, stop) {
    const exons = ivs(t.exons), cds = ivs(t.cds), mid = gy + 9;
    let g = `<g opacity="${excluded ? 0.4 : 1}" data-tip="${esc(tipText)}"><line x1="${cl(t.start)}" y1="${mid}" x2="${cl(t.end)}" y2="${mid}" stroke="${color}" stroke-width="1" ${excluded ? 'stroke-dasharray="3 2"' : ""}/>`;
    for (let px = cl(t.start) + 18; px < cl(t.end) - 10; px += 60) g += `<text x="${px}" y="${mid + 3.5}" font-size="9" fill="${color}" text-anchor="middle">${t.strand === "-" ? "‹" : "›"}</text>`;
    for (const [s, e] of exons) g += `<rect x="${cl(s)}" y="${mid - 4}" width="${Math.max(1, cl(e + 1) - cl(s))}" height="8" fill="var(--surface)" stroke="${color}" stroke-width="1.2"/>`;
    for (const [s, e] of cds) g += `<rect x="${cl(s)}" y="${mid - 7}" width="${Math.max(1, cl(e + 1) - cl(s))}" height="14" fill="${color}"/>`;
    if (stop) {
      const added = t.strand === "-" ? [stop.new_terminal_cds_start, stop.terminal_cds_start - 1] : [stop.terminal_cds_end + 1, stop.new_terminal_cds_end];
      g += `<rect x="${cl(added[0])}" y="${mid - 9}" width="${Math.max(2, cl(added[1] + 1) - cl(added[0]))}" height="18" fill="#c4271c"/>`;
    }
    g += `<text x="${4 + (markers.length ? 10 * markers.length + 2 : 0)}" y="${mid + 4}" class="mono" font-size="10">${esc(label.length > 40 ? label.slice(0, 39) + "…" : label)}</text>`;
    markers.forEach((c, i) => (g += `<rect x="${4 + i * 10}" y="${mid - 4}" width="8" height="8" fill="${c || "var(--na)"}" transform="rotate(45 ${8 + i * 10} ${mid})"/>`));
    return g + "</g>";
  }
}
async function loadGeneDetail(geneId) {
  let d;
  try { d = await api("/api/gene", { reference: S.ref, policy: S.policy, gene_id: geneId }); } catch (e) { showError($("#geneDetail"), e); return; }
  const g = d.gene, runs = visibleRuns();
  const diff = (o) => {
    if (o.cds_start_diff_bp === null || o.cds_start_diff_bp === undefined) return '<span class="na">no best CDS pair</span>';
    const p = [];
    if (o.cds_start_diff_bp) p.push(`start ${o.cds_start_diff_bp > 0 ? "+" : ""}${o.cds_start_diff_bp} bp`);
    if (o.cds_stop_diff_bp) p.push(`stop ${o.cds_stop_diff_bp > 0 ? "+" : ""}${o.cds_stop_diff_bp} bp`);
    if (o.cds_introns_missing) p.push(`${o.cds_introns_missing} reference intron(s) missing`);
    if (o.cds_introns_extra) p.push(`${o.cds_introns_extra} extra predicted intron(s)`);
    return p.length ? esc(p.join("; ")) : "identical CDS";
  };
  const rows = runs.map((r) => {
    const o = d.outcomes[r.run_id], res = result(r.run_id);
    if (!o) return `<tr>${runCell(r)}<td colspan="6" class="na">${esc(shown(res) ? "gene not in this result" : (RESULT_LABEL[res?.status] || "no result") + (res?.reason ? " — " + res.reason : ""))}</td></tr>`;
    const stop = o.query_stop_codon;
    return `<tr>${runCell(r)}<td>${outcomeChip(o, r)}<div class="cell-sub">${esc(o.classification)} / CDS ${esc(o.classification_cds)}</div></td>
      <td class="mono small">${esc(o.best_cds_ref_tx || "—")}<div class="cell-sub">${esc(o.cds_matched_id || "")} ${esc(o.best_cds_query_tx || "")}</div></td>
      <td class="small">${o.one_to_one_partner ? `<span class="status st-exact">${esc(o.one_to_one_partner)}</span>` : o.cds_exact ? '<span class="tag warn">exact, not in 1:1</span>' : '<span class="muted">—</span>'}
        ${stop ? `<div class="cell-sub">stop ${esc(stop.status)}: ${esc(stop.terminal_cds_start)}-${esc(stop.terminal_cds_end)}${stop.status === "extended" ? " → " + esc(stop.new_terminal_cds_start) + "-" + esc(stop.new_terminal_cds_end) : ""}</div>` : ""}</td>
      <td class="num">${fmtInt(o.cds_shared_bases)}<div class="cell-sub">overlap ${o.cds_overlap ?? "—"}</div></td>
      <td class="small">best pair ${esc(o.cds_intron_chain_match || "NA")} · any pair ${esc(o.cds_intron_chain_any_pair || "NA")}${o.counterpart_count >= 2 ? `<br><span class="tag warn">split</span> ${esc(o.counterpart_ids)}` : ""}</td>
      <td class="small">${diff(o)}</td></tr>`;
  }).join("");
  const iso = d.isoforms.length ? `${d.isoforms.filter((i) => i.status === "evaluated").length} evaluated of ${d.isoforms.length} isoforms: ${d.isoforms.slice(0, 8).map((i) => `${i.transcript_id} (${i.status === "evaluated" ? i.selected_by : i.exclusion_reason})`).join(", ")}${d.isoforms.length > 8 ? ", …" : ""}` : "isoform statuses not available (no evaluated-transcript export)";
  $("#geneDetail").innerHTML = `<h3>${esc(g.gene_id)} ${g.name ? "· " + esc(g.name) : ""} <span class="small muted">${esc(g.biotype || "")} · ${esc(g.chrom)}:${fmtInt(g.start)}-${fmtInt(g.end)} ${esc(g.strand)}${g.original_gene_id && g.original_gene_id !== g.gene_id ? " · file ID " + esc(g.original_gene_id) : ""}</span></h3>
    <p class="small muted">Policy ${esc(S.policy)}: ${esc(iso)}. CDS differences are the comparator's per-gene columns for the best CDS pair (positive bp = predicted CDS extends further, in the direction of transcription).</p>
    <div class="xscroll"><table class="grid compact"><tr><th>Prediction run</th><th>Outcome</th><th>Compared isoforms</th><th>1:1 exact pair</th><th class="num">Shared CDS bases</th><th>CDS intron chain</th><th>CDS differences</th></tr>${rows}</table></div>`;
}

// ---------------------------------------------------------------- provenance
function kv(obj) { return `<dl class="kv">${Object.entries(obj).map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v === null || v === undefined || v === "" ? '<span class="na">unknown</span>' : esc(typeof v === "object" ? JSON.stringify(v) : v)}</dd>`).join("")}</dl>`; }
function identityText(id, res) {
  if (!res || !shown(res)) return res ? RESULT_LABEL[res.status] : "no result";
  if (res.status === "legacy") return "pre-contract result: no implementation fingerprint (legacy import)";
  if (!id) return "unknown";
  if (id.status === "verified") return `verified implementation ${id.build_id} (historical build; ${id.scheme})`;
  if (id.status === "current_checkout") return `current code (this checkout; ${id.scheme})`;
  return `unverified implementation — ${id.reason}`;
}
function renderProvenance() {
  const m = S.meta, el = $("#provenance"), v = m.validation, ref = refMeta();
  const runIds = configRuns().map((r) => r.run_id);
  const issues = v ? v.issues.filter((i) => i.where.includes(S.ref) || runIds.some((r) => i.where.includes(r))) : [];
  let html = `<h2>Experiment</h2>${kv({ id: m.experiment.id, title: m.experiment.title, configuration: m.experiment.config_path, "configuration version": m.experiment.config_version,
    "dataset built": m.built, "store version": m.store_version, "comparator contract (this code)": `${m.contract.contract_id} ${m.contract.contract_version}, definitions ${m.contract.definition_version}`, validated: v?.validated || "not validated" })}`;
  html += `<h2>Implementation registry</h2><p class="small muted">A result's code is identified by (fingerprint scheme, fingerprint). Fingerprints of different schemes are never compared. Unknown pairs are shown as unverified implementation.</p>
    <table class="grid compact"><tr><th>Build</th><th>Commit</th><th>Scheme</th><th>Fingerprint</th><th>Files</th><th>Contract</th><th>Evidence</th></tr>${(m.registry?.builds || []).map((b) => `<tr><td>${esc(b.build_id)}</td><td class="mono">${esc((b.commit || "").slice(0, 12))}</td><td class="mono small">${esc(b.scheme)}</td><td class="mono small">${esc(b.fingerprint.slice(0, 16))}…</td><td class="num">${b.files}</td><td>${esc(b.contract_version || "")} / ${esc(b.definition_version || "")}</td><td class="small">${(b.evidence || []).map((e) => `• ${esc(e.type)}: ${esc(e.detail)}${e.sha256 ? " (" + esc(e.sha256.slice(0, 12)) + "…)" : ""}`).join("<br>")}${b.notes ? `<div class="cell-sub">${esc(b.notes)}</div>` : ""}</td></tr>`).join("")}</table>`;
  html += `<h2>Validation ${v ? (v.ok ? '<span class="tag good">no errors</span>' : '<span class="tag bad">errors</span>') : '<span class="tag warn">not run</span>'}</h2>`;
  html += issues.length ? issues.map((i) => `<div class="issue ${esc(i.level)}"><b>${esc(i.level)}</b> · ${esc(i.where)}: ${esc(i.message)}</div>`).join("") : '<p class="muted">No issues recorded for this reference and its runs.</p>';
  html += `<h2>Reference</h2>${kv({ id: ref.id, species: ref.species, assembly: ref.assembly, accession: ref.assembly_accession, source: ref.annotation_source, release: ref.annotation_release, annotation: ref.annotation, genome: ref.genome, "seqname map": ref.seqname_map, notes: ref.notes,
    "locus models": ref.models_available ? "rebuilt from the evaluated-transcript export and checked" : "unavailable: " + (ref.models_error || "") })}`;
  for (const c of ref.caveats || []) html += `<div class="issue warning">${esc(c)}</div>`;
  html += "<h2>Prediction runs</h2>";
  for (const r of orderRuns(configRuns())) {
    const ia = r.input_audit;
    html += `<h3><span class="swatch" style="background:${r.color}"></span> ${esc(runName(r))} ${isSupp(r) ? '<span class="tag info">supplementary</span>' : ""}</h3>${kv({ "run id": r.run_id, tool: r.tool, "run label": r.run_label, model: r.model, version: r.tool_version,
      variant: r.variant, "comparator options": Object.keys(r.comparator_options || {}).length ? r.comparator_options : "none (original coordinates)", "listed after": r.companion_of,
      annotation: r.annotation, "original source path": r.source?.original_path, "source sha256": r.source?.sha256, "recorded in": r.source?.recorded_in, "legacy ids": Object.keys(r.legacy_ids || {}).length ? r.legacy_ids : "",
      "locus models": m.model_errors?.["run:" + r.run_id] ? "unavailable: " + m.model_errors["run:" + r.run_id] : `rebuilt (${fmtInt(m.model_components?.["run:" + r.run_id]?.transcripts)} transcripts; ${fmtInt(m.model_components?.["run:" + r.run_id]?.transformed)} with an audited transformation)`,
      "sequences with prediction records": ia?.sequences_with_records, "processed sequences recorded": ia ? (ia.processed_sequences_record?.sequences ? `${ia.processed_sequences_record.sequences.length ?? ia.processed_sequences_record.sequences} (${ia.processed_sequences_record.source})` : `not recorded (${ia.processed_sequences_record?.source || "unknown"})`) : "",
      "reference PC genes on sequences without records": ia?.reference_pc_genes_on_sequences_without_records ? `${fmtInt(ia.reference_pc_genes_on_sequences_without_records.genes)} on ${fmtInt(ia.reference_pc_genes_on_sequences_without_records.sequences)} sequences (counted as missed; not established as unprocessed)` : "" })}`;
    for (const c of r.caveats || []) html += `<div class="issue warning">${esc(c)}</div>`;
    html += `<table class="grid compact"><tr><th>Policy</th><th>Status</th><th>Implementation</th><th>Output checksums</th><th>1:1 pairs file</th><th>Manifest provenance</th><th>Output directory</th></tr>${m.policies.map((p) => {
      const res = result(r.run_id, p.id) || {};
      const prov = res.provenance || {};
      return `<tr><td>${esc(p.id)}</td><td>${esc(RESULT_LABEL[res.status] || res.status || "no result")}${res.reason && res.status !== "complete" ? `<div class="cell-sub">${esc(res.reason)}</div>` : ""}</td>
        <td class="small">${esc(identityText(res.identity, res))}${res.identity?.fingerprint ? `<div class="cell-sub mono">${esc(res.identity.fingerprint.slice(0, 16))}…</div>` : ""}</td>
        <td class="small">${res.outputs_check ? `${esc(res.outputs_check.status)} (${res.outputs_check.files} files)` : ""}</td>
        <td class="small">${res.pairs_check ? `${esc(res.pairs_check.status)}${res.pairs_check.rows !== undefined ? ` (${fmtInt(res.pairs_check.rows)} rows, TP ${fmtInt(res.pairs_check.true_positives)})` : ""}` : ""}</td>
        <td class="small">${prov.prediction_id ? `${esc(prov.prediction_id)} · ${esc(prov.predictor || "")}${prov.run_label ? " · " + esc(prov.run_label) : ""} · ${esc(prov.variant || "")}` : ""}</td>
        <td class="mono small">${esc(res.output_dir || "")}</td></tr>`;
    }).join("")}</table>`;
  }
  const imp = m.import_report;
  if (imp) {
    html += `<h2>Import</h2>${kv({ directories: (imp.directories || []).map((d) => `${d.kind}: ${d.path}`).join("; "), imported: imp.imported, "results matched": imp.matched.length, "pilots (kept separate)": imp.pilots.length, unmatched: imp.unmatched.length, conflicts: imp.conflicts.length })}`;
    if (imp.pilots.length) html += `<details class="panel"><summary>Pilots — not headline results, never shown in the overview</summary><table class="grid compact"><tr><th>Directory</th><th>Run</th><th>Policy</th><th>Region</th><th>Implementation</th><th>Reason</th></tr>${imp.pilots.map((p) => `<tr><td class="mono small">${esc(p.dir)}</td><td>${esc(p.run_id)}</td><td>${esc(p.policy_id || "")}</td><td>${esc((p.options_region || []).join(","))}</td><td class="small">${esc(p.implementation || "")}</td><td class="small">${esc(p.reason)}</td></tr>`).join("")}</table></details>`;
    if (imp.unmatched.length) html += `<details class="panel"><summary>Unmatched directories</summary><pre class="mono">${esc(JSON.stringify(imp.unmatched, null, 1))}</pre></details>`;
  }
  el.innerHTML = html;
}

// ---------------------------------------------------------------- dictionary
function renderDictionary() {
  const d = S.meta.metric_dictionary, defs = S.meta.metric_definitions ? Object.values(S.meta.metric_definitions)[0] || {} : {};
  $("#dictionary").innerHTML = `<h2>Metric dictionary</h2><ul>${d.notes.map((n) => `<li class="small">${esc(n)}</li>`).join("")}</ul>
    <table class="grid compact"><tr><th>Metric</th><th>Group</th><th>Numerator (manifest)</th><th>Denominator</th><th>Unit</th><th>Explanation and limitations</th></tr>
    ${d.metrics.map((m) => { const def = defs[m.id] || {}; return `<tr><td><b>${esc(m.label)}</b><div class="cell-sub mono">${esc(m.id)}</div></td><td>${esc(m.group)}</td><td class="small">${esc(def.numerator || "")}</td><td class="small mono">${esc(def.denominator_id || "")}</td><td>${esc(def.unit || "")}</td>
      <td class="small">${esc(m.definition)}${(m.limitations || []).map((l) => "<br>• " + esc(l)).join("")}</td></tr>`; }).join("")}</table>`;
}

window.addEventListener("hashchange", () => {
  const p = new URLSearchParams(location.hash.slice(1));
  if (p.get("exp") && p.get("exp") !== S.exp) { S.exp = p.get("exp"); $("#exp").value = S.exp; loadExperiment(p); return; }
  let changed = false;
  if (p.get("policy") && p.get("policy") !== S.policy && S.meta.policies.some((x) => x.id === p.get("policy"))) { S.policy = p.get("policy"); changed = true; }
  if (p.get("scope") && p.get("scope") !== S.scope && scopes().some((x) => x.id === p.get("scope") && x.available)) { S.scope = p.get("scope"); changed = true; }
  if (changed) { S.genesLoaded = S.predsLoaded = false; renderAll(); }
});

init().catch((e) => { document.body.insertAdjacentHTML("afterbegin", `<div class="issue error">${esc(e.message)}</div>`); });
