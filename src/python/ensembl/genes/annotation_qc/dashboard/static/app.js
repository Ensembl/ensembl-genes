"use strict";
// Annotation QC dashboard (no external dependencies).
const $ = (s, el = document) => el.querySelector(s);
const $$ = (s, el = document) => [...el.querySelectorAll(s)];
const esc = (v) => String(v ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const fmtInt = (n) => (n === null || n === undefined || Number.isNaN(n) ? "—" : Math.round(n).toLocaleString());
const fmtPct = (v, d = 1) => (v === null || v === undefined ? "—" : (100 * v).toFixed(d) + "%");

const STATUS = {
  exact: ["exact", "coordinate-exact CDS"],
  terminal_diff: ["ends ≠", "same CDS intron chain; start/stop differs"],
  diff_chain: ["splice ≠", "≥ 0.8 reciprocal CDS overlap, different CDS intron chain"],
  partial_lt08: ["< 0.8", "shares CDS, < 0.8 reciprocal overlap"],
  span_only: ["span", "same-strand gene-span partner without CDS overlap"],
  strand_mismatch: ["strand", "opposite-strand overlap only"],
  missed: ["missed", "no same-strand partner"],
};
const RESULT_LABEL = { complete: "complete", failed: "failed", not_run: "not run", stale: "stale", unavailable: "unavailable" };
const HEADLINE = ["locus_recovery_cds_overlap", "cds_exact_ref", "cds_exact_query", "exact_recall_one_to_one", "exact_precision_one_to_one",
  "exact_f1_one_to_one", "cds_intron_chain_best_pair", "cds_intron_chain_any_pair", "ref_cds_introns_recovered", "query_cds_introns_supported",
  "missed_reference_genes", "novel_predictions", "splits", "merges"];
const STRUCTURE = ["exact_tp_one_to_one", "exact_unmatched_reference_genes", "exact_unmatched_query_genes", "exact_gene_pairs",
  "exact_components_with_several_edges", "cds_exact_ref_multi", "cds_exact_ref_single", "partial_cds", "partial_same_chain_terminal_diff",
  "partial_ge08_diff_chain", "partial_lt08", "cds_structural_exact_match", "gene_span_locus_detection", "gene_span_only_no_exon_overlap",
  "detected_without_cds_overlap", "query_matched_without_cds_overlap", "exact_pairs_query_gene_reused"];
const STOP = ["stop_codon_transcripts_with_cds", "stop_codon_eligible", "stop_codon_extended", "stop_codon_unchanged", "stop_codon_skipped"];
const HISTORICAL = ["historical_exact_tp_best_pair", "historical_exact_recall_best_pair", "historical_exact_precision_best_pair", "historical_exact_f1_best_pair"];
const METRIC_STATUS = { not_available: "n/a", undefined: "undefined", not_applicable: "—", inconsistent: "inconsistent" };
const LOCI = ["missed_reference_genes", "ref_genes_without_cds_overlap", "novel_predictions", "splits", "merges", "one_to_one_ref_genes",
  "strand_mismatch_ref", "strand_mismatch_query"];
const NOVEL = ["novel_Overlaps_reference_pseudogene", "novel_Overlaps_other_non_protein_coding_reference",
  "novel_Overlaps_unevaluated_protein_coding_reference", "novel_Overlaps_evaluated_reference_opposite_strand", "novel_Novel_no_reference_overlap"];
const NOVEL_COLORS = ["#184f95", "#2a78d6", "#6da7ec", "#9ec5f4", "#cde2fb"];
const UTR = ["exon_coordinate_exact_ref", "structural_exact_match_with_utr"];
const SCALE = ["reference_coding_genes", "reference_eligible_transcripts", "predicted_coding_genes", "predicted_transcripts"];

const WHOLE = "whole_assembly";
const S = { exp: null, meta: null, ref: null, policy: null, scope: WHOLE, view: "coverage", runs: new Set(), ov: null, dict: {}, genePage: 0, predPage: 0,
  locus: null, selectedGene: null, pageSize: 100, plotPolicy: null, posCoord: "cds_span", posMin: 0, lenMeasure: "exact", splice: {} };

async function api(path, params = {}) {
  const url = new URL(path, location.origin);
  for (const [k, v] of Object.entries({ exp: S.exp, ...params })) if (v !== undefined && v !== null && v !== "") url.searchParams.set(k, v);
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || res.statusText);
  return data;
}
function apiUrl(path, params) {
  const url = new URL(path, location.origin);
  for (const [k, v] of Object.entries({ exp: S.exp, ...params })) if (v !== undefined && v !== null && v !== "") url.searchParams.set(k, v);
  return url.toString();
}
function showError(where, err) {
  where.innerHTML = `<div class="issue error">${esc(err.message || err)}</div>`;
}

// ---------------------------------------------------------------- tooltip
const tip = $("#tooltip");
document.addEventListener("mousemove", (e) => {
  const t = e.target.closest("[data-tip]");
  if (!t) { tip.hidden = true; return; }
  tip.textContent = t.dataset.tip;
  tip.hidden = false;
  const x = Math.min(e.clientX + 14, window.innerWidth - 380);
  tip.style.left = x + "px"; tip.style.top = (e.clientY + 14) + "px";
});

// ---------------------------------------------------------------- setup
async function init() {
  const exps = await api("/api/experiments");
  $("#exp").innerHTML = exps.map((e) => `<option value="${esc(e.id)}">${esc(e.title || e.id)}${e.synthetic ? " [SYNTHETIC]" : ""}</option>`).join("");
  const params = new URLSearchParams(location.hash.slice(1));
  S.exp = params.get("exp") && exps.some((e) => e.id === params.get("exp")) ? params.get("exp") : exps[0].id;
  $("#exp").value = S.exp;
  $("#exp").onchange = () => { S.exp = $("#exp").value; loadExperiment(); };
  $$(".tabs button").forEach((b) => (b.onclick = () => selectTab(b.dataset.tab)));
  $$(".subnav button").forEach((b) => (b.onclick = () => selectView(b.dataset.view)));
  await loadExperiment(params);
}
function selectView(name) {
  S.view = name;
  $$(".subnav button").forEach((b) => b.classList.toggle("active", b.dataset.view === name));
  $$(".view").forEach((v) => v.classList.toggle("active", v.id === "view-" + name));
  if (name === "structure") renderPlots();
  if (name === "detail") { renderPolicyChart(); renderPartialChart(); renderNovelChart(); }
  writeHash();
}
function selectTab(name) {
  $$(".tabs button").forEach((b) => b.classList.toggle("active", b.dataset.tab === name));
  $$(".tab").forEach((t) => t.classList.toggle("active", t.id === "tab-" + name));
  if (name === "genes" && !S.genesLoaded) loadGenes();
  if (name === "predictions" && !S.predsLoaded) loadPredictions();
}
async function loadExperiment(params = new URLSearchParams()) {
  S.runs.clear(); // run ids belong to one experiment
  S.meta = await api("/api/meta");
  S.dict = Object.fromEntries((S.meta.metric_dictionary?.metrics || []).map((m) => [m.id, m]));
  $("#synthetic").hidden = !S.meta.experiment.synthetic;
  const refs = S.meta.references;
  $("#ref").innerHTML = refs.map((r) => `<option value="${esc(r.id)}">${esc(r.species)} · ${esc(r.assembly)} · ${esc(r.id)}</option>`).join("");
  S.ref = params.get("ref") && refs.some((r) => r.id === params.get("ref")) ? params.get("ref") : refs[0].id;
  $("#ref").value = S.ref;
  $("#ref").onchange = () => { S.ref = $("#ref").value; S.runs.clear(); loadReference(); };
  S.policy = params.get("policy") || S.meta.policies[0].id;
  S.scope = params.get("scope") || WHOLE;
  S.plotPolicy = params.get("plot") || S.plotPolicy;
  if (["cds_span", "coding"].includes(params.get("pos"))) S.posCoord = params.get("pos");
  if (["0", "60"].includes(params.get("posmin"))) S.posMin = Number(params.get("posmin"));
  if (["exact", "donor", "acceptor"].includes(params.get("len"))) S.lenMeasure = params.get("len");
  S.splice = {};
  renderDictionary(); renderContext();
  await loadReference();
  if (params.get("view")) selectView(params.get("view"));
}
function configRuns() { return S.meta.runs.filter((r) => r.reference === S.ref); }
// Display order: pinned runs first, the others by the configured metric (descending) for the
// selected policy and scope; a supplementary companion follows its run and is never ranked.
// Without a display configuration the configured order is kept.
function orderRuns(runs) {
  const d = S.meta.display || {};
  const configured = Object.keys(d).length > 0;
  const byOrder = (a, b) => a.config_order - b.config_order;
  const companion = (r) => r.companion_of && runs.some((x) => x.run_id === r.companion_of);
  const leaders = runs.filter((r) => !companion(r));
  const pinned = (d.pinned_runs || []).map((id) => leaders.find((r) => r.run_id === id)).filter(Boolean);
  const value = (r) => { const m = S.ov ? scopedMetric(r.run_id, S.policy, d.order_metric || "exact_recall_one_to_one") : null; return m && m.status === "ok" && m.value !== null ? m.value : -1; };
  const rest = leaders.filter((r) => !pinned.includes(r)).sort((a, b) => !configured ? byOrder(a, b)
    : (a.role === "supplementary") - (b.role === "supplementary") || value(b) - value(a) || byOrder(a, b));
  const out = [];
  for (const r of [...pinned, ...rest]) { out.push(r); runs.filter((c) => companion(c) && c.companion_of === r.run_id).sort(byOrder).forEach((c) => out.push(c)); }
  return out;
}
function refRuns() { return orderRuns(configRuns()); }
function refMeta() { return S.meta.references.find((r) => r.id === S.ref); }
function scopes() { return refMeta().scopes || [{ id: WHOLE, label: "All compared sequences", available: true }]; }
function scopeInfo() { return scopes().find((x) => x.id === S.scope) || scopes()[0]; }
function scopedMetric(run, policy, id) {
  if (S.scope === WHOLE) return metric(run, policy, id);
  return (S.ov.scope_metrics || []).find((m) => m.scope_id === S.scope && m.run_id === run && m.policy_id === policy && m.metric_id === id);
}
function runDisplay(r) {
  const same = configRuns().filter((x) => runName(x) === runName(r)).length > 1;
  return same ? `${runName(r)} (${r.run_id})` : runName(r);
}
function runById(id) { return S.meta.runs.find((r) => r.run_id === id); }
function selectedRuns() { return refRuns().filter((r) => S.runs.has(r.run_id)); }
function result(run, policy) { return S.ov.results.find((r) => r.run_id === run && r.policy_id === policy); }
function metric(run, policy, id) { return S.ov.metrics.find((m) => m.run_id === run && m.policy_id === policy && m.metric_id === id); }
function runName(r) { return r.label || r.tool; }
function harmonised(r) { return r.query_coordinates === "stop_harmonised"; }
function runTags(r) {
  return `${r.role === "supplementary" ? '<span class="tag info" data-tip="supplementary result: not a headline comparison">supplementary</span>' : ""}${harmonised(r)
    ? `<span class="tag warn" data-tip="${esc("query compared after --add-stop-codon (genetic code " + (r.comparator_options?.stop_codon_genetic_code ?? "?") + "): terminal CDS extended by an adjacent in-frame stop where the frame is established from the sequence. Locus tracks show these evaluated coordinates.")}">stop-harmonised</span>` : ""}`;
}
function runLabel(r, html = true) {
  const detail = [r.model ? `model ${r.model}` : "model unknown", r.tool_version ? `v${r.tool_version}` : "version unknown",
    harmonised(r) ? "query coordinates: stop-harmonised (--add-stop-codon)" : "query coordinates: original"].join(", ");
  if (!html) return runDisplay(r);
  return `<span class="swatch" style="background:${r.color}"></span><span data-tip="${esc("run id " + r.run_id + "\n" + detail)}"><b>${esc(runDisplay(r))}</b></span> ${runTags(r)}`;
}

async function loadReference() {
  S.ov = await api("/api/overview", { reference: S.ref });
  S.ov.scope_metrics = S.ov.scope_metrics || []; S.ov.analyses = S.ov.analyses || [];
  if (!scopes().some((x) => x.id === S.scope && x.available)) S.scope = WHOLE;
  if (!S.runs.size) refRuns().forEach((r) => S.runs.add(r.run_id));
  renderPolicySelector(); renderScopeSelector(); renderRuns(); renderOverview(); renderProvenance();
  S.genesLoaded = S.predsLoaded = false;
  $("#locusPanel").hidden = true;
  setupGeneFilters(); setupPredFilters();
  const active = $(".tabs button.active").dataset.tab;
  if (active === "genes") loadGenes();
  if (active === "predictions") loadPredictions();
  writeHash();
}
function writeHash() {
  history.replaceState(null, "", `#exp=${encodeURIComponent(S.exp)}&ref=${encodeURIComponent(S.ref)}&policy=${encodeURIComponent(S.policy)}`
    + `&scope=${encodeURIComponent(S.scope)}&view=${encodeURIComponent(S.view)}`
    + (S.view === "structure" ? `&plot=${encodeURIComponent(S.plotPolicy || "")}&pos=${S.posCoord}&posmin=${S.posMin}&len=${S.lenMeasure}` : ""));
}
function policyAvailability(pid) {
  const rs = refRuns().map((r) => result(r.run_id, pid)).filter(Boolean);
  if (!rs.length) return { ok: false, why: "no results" };
  if (rs.every((r) => r.status === "unavailable")) return { ok: false, why: rs[0].reason };
  return { ok: true };
}
function renderPolicySelector() {
  const box = $("#policy");
  box.innerHTML = S.meta.policies.map((p) => {
    const av = policyAvailability(p.id);
    return `<button role="radio" data-p="${esc(p.id)}" class="${p.id === S.policy ? "active" : ""} ${av.ok ? "" : "unavailable"}"
      data-tip="${esc(p.label + "\n" + (p.description || "") + (p.role === "supplementary" ? "\nSupplementary view" : "") + (av.ok ? "" : "\nUnavailable: " + av.why))}">${esc(p.id)} · ${esc(p.label)}${p.role === "supplementary" ? " (suppl.)" : ""}</button>`;
  }).join("");
  $$("button", box).forEach((b) => (b.onclick = () => { S.policy = b.dataset.p; renderPolicySelector(); renderRuns(); renderOverview(); S.genesLoaded = S.predsLoaded = false;
    const t = $(".tabs button.active").dataset.tab; if (t === "genes") loadGenes(); if (t === "predictions") loadPredictions(); if (S.locus) drawLocus(); writeHash(); }));
  const p = S.meta.policies.find((x) => x.id === S.policy);
  const av = policyAvailability(S.policy);
  $("#policyNote").innerHTML = `${esc(p ? p.description : "")}${av.ok ? "" : ` <span class="tag warn">unavailable for this reference: ${esc(av.why)}</span>`}`;
}
function renderScopeSelector() {
  const sel = $("#scope");
  sel.innerHTML = scopes().map((x) => `<option value="${esc(x.id)}" ${x.available ? "" : "disabled"}>${esc(x.label)}${x.role === "supplementary" ? " (suppl.)" : ""}${x.available ? "" : " — unavailable"}</option>`).join("");
  sel.value = S.scope;
  sel.disabled = scopes().length < 2;
  sel.title = scopeInfo().description || scopeInfo().reason || "";
  sel.onchange = () => { S.scope = sel.value; renderRuns(); renderOverview(); writeHash(); };
}
function renderRuns() {
  $("#runs").innerHTML = refRuns().map((r) => `<label class="chip"><input type="checkbox" data-run="${esc(r.run_id)}" ${S.runs.has(r.run_id) ? "checked" : ""}>${runLabel(r)}</label>`).join(" ")
    || '<span class="na">no runs configured for this reference</span>';
  $$("#runs input").forEach((c) => (c.onchange = () => { c.checked ? S.runs.add(c.dataset.run) : S.runs.delete(c.dataset.run); renderOverview();
    if (S.genesLoaded) renderGeneTable(S.lastGenes); if (S.locus) drawLocus(); }));
}

// ---------------------------------------------------------------- overview
function metricCell(runId, policy, id) {
  const res = result(runId, policy);
  if (!res) return `<td class="num"><span class="na">no result</span></td>`;
  if (res.status !== "complete") return `<td class="num"><span class="na" data-tip="${esc(res.reason || "")}">${esc(RESULT_LABEL[res.status] || res.status)}</span></td>`;
  const m = metric(runId, policy, id);
  if (!m) return `<td class="num"><span class="na" data-tip="metric not in this dataset">—</span></td>`;
  const def = S.dict[id] || {};
  const src = m.source && def.source && m.source !== def.source ? `<div class="cell-sub"><span class="tag">${esc(m.source)}</span></div>` : "";
  if (m.status !== "ok" && m.status !== "inconsistent")
    return `<td class="num"><span class="na" data-tip="${esc(m.status + (m.note ? ": " + m.note : ""))}">${esc(METRIC_STATUS[m.status] || m.status)}</span></td>`;
  const warn = m.status === "inconsistent" ? `<div class="cell-sub"><span class="tag bad" data-tip="${esc(m.note || "")}">inconsistent</span></div>` : "";
  const tip = m.note ? ` data-tip="${esc(m.note)}"` : "";
  if (def.unit === "rate") {
    if (m.value === null || m.value === undefined) return `<td class="num"><span class="na" data-tip="${esc(m.note || "no value")}">n/a</span></td>`;
    const parts = id === "exact_f1_one_to_one" ? `2TP ${fmtInt(m.numerator)} / R+Q ${fmtInt(m.denominator)}` : `${fmtInt(m.numerator)}${m.denominator ? " / " + fmtInt(m.denominator) : ""}`;
    return `<td class="num"${tip}><div class="cell-main">${m.value.toFixed(4)}</div><div class="cell-sub">${parts}</div>${src}${warn}</td>`;
  }
  if (m.denominator === null || m.denominator === undefined) return `<td class="num"${tip}><div class="cell-main">${fmtInt(m.numerator)}</div>${src}${warn}</td>`;
  return `<td class="num"${tip}><div class="cell-main">${fmtPct(m.numerator / m.denominator)}</div><div class="cell-sub">${fmtInt(m.numerator)} / ${fmtInt(m.denominator)}</div>${src}${warn}</td>`;
}
function metricTable(ids, el) {
  const runs = selectedRuns();
  if (!runs.length) { el.innerHTML = '<p class="na">Select at least one prediction run.</p>'; return; }
  const head = `<tr><th>Metric</th><th>Denominator</th><th>Source</th>${runs.map((r) => `<th>${runLabel(r)}</th>`).join("")}</tr>`;
  const body = ids.map((id) => {
    const d = S.dict[id] || { label: id };
    const t = [d.definition, "Numerator: " + (d.numerator || ""), "Denominator: " + (d.denominator || ""), ...(d.limitations || []).map((l) => "Limitation: " + l)].filter(Boolean).join("\n");
    return `<tr><td><span class="info" data-tip="${esc(t)}">${esc(d.label || id)}</span></td><td class="small muted">${esc(d.denominator || "")}</td>
      <td><span class="tag ${d.source === "independent" ? "info" : ""}">${esc(d.source || "")}</span></td>${runs.map((r) => metricCell(r.run_id, S.policy, id)).join("")}</tr>`;
  }).join("");
  el.innerHTML = `<table class="grid">${head}${body}</table>
    <div class="small"><button class="btn" data-dl>Download this table (TSV)</button></div>`;
  $("[data-dl]", el).onclick = () => downloadTSV(`${S.ref}_${S.policy}_metrics.tsv`, tableRows(ids, runs));
}
function tableRows(ids, runs) {
  const rows = [["metric_id", "label", "source", "run_id", "tool", "policy", "result_status", "metric_status", "numerator", "denominator", "value", "note"]];
  for (const id of ids) for (const r of runs) {
    const res = result(r.run_id, S.policy); const m = metric(r.run_id, S.policy, id);
    rows.push([id, S.dict[id]?.label || id, m?.source || S.dict[id]?.source || "", r.run_id, runName(r), S.policy, res?.status || "no_result", m?.status || "",
      m?.numerator ?? "", m?.denominator ?? "", m?.value ?? "", m?.note || ""]);
  }
  return rows;
}
function downloadTSV(name, rows) {
  const blob = new Blob([rows.map((r) => r.map((v) => String(v).replace(/\t|\n/g, " ")).join("\t")).join("\n") + "\n"], { type: "text/tab-separated-values" });
  const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
function renderOverview() {
  renderContextBar(); renderSummary(); renderStructureSummary(); renderDetail();
  if (S.view === "structure") renderPlots();
}
function renderDetail() {
  const ref = S.meta.references.find((r) => r.id === S.ref);
  $("#detailScopeNote").hidden = S.scope === WHOLE;
  $("#detailScopeNote").textContent = `Detailed metrics are the comparator's whole-assembly values; the scope "${scopeInfo().label}" applies to the coverage and structural views.`;
  const prep = ref.preparation || {};
  const R = refRuns().map((r) => metric(r.run_id, S.policy, "reference_coding_genes")).find((m) => m && m.status === "ok");
  const scope = prep.scope ? (prep.scope.type === "all" ? "all sequences" : `${(prep.scope.sequences || []).length} sequences (${prep.scope.type})`) : "unknown (not prepared)";
  $("#legacy").hidden = !S.meta.legacy_dataset;
  $("#legacy").textContent = S.meta.legacy_note || "";
  $("#scopeLine").innerHTML = `<b>${esc(ref.species)}</b> · ${esc(ref.assembly)} ${ref.assembly_accession ? "(" + esc(ref.assembly_accession) + ")" : ""} ·
    ${esc(ref.annotation_source || "")} ${esc(ref.annotation_release || "")}<br>
    Scope: ${esc(scope)} · Evaluation: ${esc(S.meta.evaluation.evaluation_mode)}, transcript biotypes ${esc((S.meta.evaluation.reference_transcript_biotypes || []).join(",") || "any")} ·
    R = ${R ? fmtInt(R.numerator) : "—"} reference coding genes under policy ${esc(S.policy)}${ref.models_available ? "" : ' · <span class="tag warn">reference models unavailable: locus browser disabled</span>'}`;
  $("#statusStrip").innerHTML = `<div class="strip">${refRuns().map((r) => {
    const res = result(r.run_id, S.policy) || { status: "no result" };
    const ind = res.independent?.status;
    const cls = res.status === "complete" ? "good" : res.status === "unavailable" ? "info" : "bad";
    return `<div class="item">${runLabel(r)}<br><span class="tag ${cls}">${esc(RESULT_LABEL[res.status] || res.status)}</span>
      ${res.source ? `<span class="tag">${esc(res.source)}</span>` : ""} ${ind ? `<span class="tag ${ind === "verified" ? "good" : "info"}" data-tip="independent analysis: ${esc(ind)}">indep. ${esc(ind)}</span>` : ""}
      ${(r.caveats || []).length ? `<span class="tag warn" data-tip="${esc(r.caveats.join("\n"))}">${r.caveats.length} caveat${r.caveats.length > 1 ? "s" : ""}</span>` : ""}
      ${res.status !== "complete" && res.reason ? `<div class="small muted">${esc(res.reason)}</div>` : ""}</div>`;
  }).join("")}</div>`;
  metricTable(HEADLINE, $("#headlineTable"));
  metricTable(STRUCTURE, $("#structureTable"));
  metricTable(LOCI, $("#lociTable"));
  metricTable(NOVEL, $("#novelTable"));
  metricTable(UTR, $("#utrTable"));
  metricTable(SCALE, $("#scaleTable"));
  metricTable(STOP, $("#stopTable"));
  metricTable(HISTORICAL, $("#historicalTable"));
  renderDownloads();
  const sel = $("#policyMetric");
  const rateIds = HEADLINE.concat(STRUCTURE).filter((id) => S.dict[id] && S.dict[id].denominator && S.dict[id].denominator !== "—");
  if (!sel.options.length) { sel.innerHTML = rateIds.map((id) => `<option value="${id}">${esc(S.dict[id].label)}</option>`).join(""); sel.value = "exact_f1_one_to_one"; sel.onchange = renderPolicyChart; }
  if (S.view === "detail") { renderPolicyChart(); renderPartialChart(); renderNovelChart(); }
  $("#dlMetrics").onclick = () => (location.href = apiUrl("/api/overview", { reference: S.ref, format: "csv" }));
}

function renderDownloads() {
  const runs = selectedRuns().filter((r) => result(r.run_id, S.policy)?.status === "complete");
  const el = $("#downloads");
  if (!runs.length || S.meta.legacy_dataset) { el.innerHTML = S.meta.legacy_dataset ? '<span class="na">one-to-one pairs are not stored in this (store version 1) dataset</span>' : ""; return; }
  el.innerHTML = `<span class="small muted">One-to-one exact-CDS pairs (policy ${esc(S.policy)}, TSV):</span> ${runs.map((r) => {
    const chk = result(r.run_id, S.policy)?.pairs_check || {};
    return `<a class="btn" href="${esc(apiUrl("/api/pairs", { reference: S.ref, policy: S.policy, run: r.run_id, format: "csv" }))}" data-tip="${esc(`${chk.rows ?? "?"} pairs exported; TP ${chk.true_positives ?? "?"}; check: ${chk.status || "?"}`)}">${esc(runName(r))}</a>`;
  }).join(" ")}${runs.filter(harmonised).map((r) => ` <a class="btn" href="${esc(apiUrl("/api/stop_audit", { run: r.run_id, format: "csv" }))}">stop-codon audit: ${esc(runName(r))}</a>`).join("")}`;
}

// ---------------------------------------------------------------- summary views
const SELECTION_TEXT = { all: "any eligible isoform", longest_cds: "the longest CDS per gene", canonical: "the canonical transcript" };
const fmtF1 = (v) => (v === null || v === undefined ? "—" : v.toFixed(3));
function policyInfo(id = S.policy) { return S.meta.policies.find((p) => p.id === id); }
function policySentence(p) {
  if (!p) return "";
  if (p.reference_transcript_selection === "all")
    return `View ${p.id} (${p.label}): a reference gene counts when any of its eligible isoforms is reproduced — locus / any-isoform recovery.`;
  return `View ${p.id} (${p.label}): scores ${SELECTION_TEXT[p.reference_transcript_selection] || p.reference_transcript_selection} on the reference side and ${SELECTION_TEXT[p.query_transcript_selection] || p.query_transcript_selection} on the predicted side, one model per gene; a correct model of another isoform counts as a miss.`;
}
function scopeR() {
  const m = refRuns().map((r) => scopedMetric(r.run_id, S.policy, "reference_coding_genes")).find((x) => x && x.status === "ok");
  return m ? m.numerator : null;
}
function renderContextBar() {
  const p = policyInfo(), sc = scopeInfo(), R = scopeR();
  const stop = refRuns().filter(harmonised);
  const other = S.meta.policies.filter((x) => x.id !== S.policy && x.reference_transcript_selection !== "all" && x.role !== "supplementary")[0];
  $("#contextBar").innerHTML = `<span class="ctx"><b>Transcript policy</b> View ${esc(p.id)} · ${esc(p.label)}${p.role === "supplementary" ? ' <span class="tag info">supplementary</span>' : ""}</span>
    <span class="ctx" data-tip="${esc(sc.description || "")}"><b>Sequence scope</b> ${esc(sc.label)}${sc.sequence_count ? ` (${fmtInt(sc.sequence_count)} sequences)` : ""}${sc.role === "supplementary" ? ' <span class="tag info">supplementary</span>' : ""} · R = ${R === null ? "—" : fmtInt(R)} reference genes</span>
    <span class="ctx"><b>Stop codons</b> ${stop.length ? `original coordinates; ${stop.map((r) => esc(runDisplay(r))).join(", ")} adds the adjacent in-frame stop (supplementary row)` : "original coordinates for every run"}</span>`;
  $("#policyExplain").textContent = policySentence(p) + (p.reference_transcript_selection === "all" && other ? ` ${policySentence(other).replace(/^View/, "View")}` : "");
}
function statusOf(run) {
  const res = result(run.run_id, S.policy);
  if (!res) return { ok: false, text: "no result" };
  if (res.status !== "complete") return { ok: false, text: RESULT_LABEL[res.status] || res.status, tip: res.reason };
  return { ok: true };
}
function naCell(text, tip, cls = "num") { return `<td class="${cls}"><span class="na" data-tip="${esc(tip || "")}">${esc(text)}</span></td>`; }
function valueCell(run, id, fmt, tipFn) {
  const st = statusOf(run);
  if (!st.ok) return naCell(st.text, st.tip);
  const m = scopedMetric(run.run_id, S.policy, id);
  if (!m) return naCell("not available", "metric not in this dataset");
  if (m.status !== "ok" && m.status !== "inconsistent")
    return naCell(m.status === "not_applicable" ? "—" : m.status === "undefined" ? "undefined" : "not available", `${m.status}${m.note ? ": " + m.note : ""}`);
  // exact ratio when counts are stored (avoids re-rounding a 4-dp value)
  const v = m.denominator && m.numerator !== null ? m.numerator / m.denominator : m.value;
  if (v === null || v === undefined) return naCell("not available", m.note || "no value");
  const shown = fmt === "f1" ? (v > 0 && v < 0.001 ? "< 0.001" : fmtF1(v)) : v > 0 && v < 0.001 ? "< 0.1%" : fmtPct(v);
  const warn = m.status === "inconsistent" ? ` <span class="tag bad" data-tip="${esc(m.note || "")}">inconsistent</span>` : "";
  const bar = `<div class="minibar"><span style="width:${(100 * Math.min(1, Math.max(0, v))).toFixed(1)}%;background:${run.color}"></span></div>`;
  return `<td class="num" data-tip="${esc(tipFn(m, v))}"><div class="cell-main">${shown}${warn}</div>${bar}</td>`;
}
function runCell(r) {
  const comp = r.companion_of ? '<span class="muted">↳ </span>' : "";
  const cav = (r.caveats || []).length ? ` <span class="caveat" data-tip="${esc("Caveats:\n• " + r.caveats.join("\n• "))}">ⓘ</span>` : "";
  return `<td class="runcell">${comp}${runLabel(r)}${cav}</td>`;
}
const SUMMARY_COLUMNS = [
  { id: "locus_recovery_cds_overlap", label: "Loci recovered", sub: "same-strand CDS overlap · % of R", fmt: "pct",
    tip: (m, v) => `${fmtInt(m.numerator)} of ${fmtInt(m.denominator)} reference genes share at least one CDS base with a same-strand prediction (${fmtPct(v, 2)})` },
  { id: "exact_recall_one_to_one", label: "Exact CDS recall", sub: "% of R in a one-to-one exact pair", fmt: "pct",
    tip: (m, v) => `TP ${fmtInt(m.numerator)} one-to-one coordinate-exact pairs / R ${fmtInt(m.denominator)} reference genes = ${fmtPct(v, 2)}` },
  { id: "exact_precision_one_to_one", label: "Exact CDS precision", sub: "% of Q in a one-to-one exact pair", fmt: "pct",
    tip: (m, v) => `TP ${fmtInt(m.numerator)} / Q ${fmtInt(m.denominator)} predicted genes = ${fmtPct(v, 2)}` },
  { id: "exact_f1_one_to_one", label: "Exact CDS F1", sub: "one-to-one · 0–1 scale", fmt: "f1",
    tip: (m, v) => `2TP ${fmtInt(m.numerator)} / (R + Q) ${fmtInt(m.denominator)} = ${v.toFixed(4)}` },
];
function columnHead(c, def) {
  const d = S.dict[def] || {};
  const t = [d.definition, d.numerator ? "Numerator: " + d.numerator : "", d.denominator ? "Denominator: " + d.denominator : ""].filter(Boolean).join("\n");
  return `<th class="num"><span class="info" data-tip="${esc(t)}">${esc(c.label)}</span><div class="cell-sub">${esc(c.sub)}</div></th>`;
}
function renderSummary() {
  const runs = selectedRuns();
  const bad = runs.filter((r) => !statusOf(r).ok);
  $("#resultStatus").innerHTML = bad.length ? `<div class="issue warning">${bad.map((r) => `${esc(runDisplay(r))}: ${esc(statusOf(r).text)}`).join(" · ")}</div>` : "";
  if (!runs.length) { $("#summaryTable").innerHTML = '<p class="na">Select at least one prediction run.</p>'; $("#summaryCounts").innerHTML = ""; return; }
  const head = `<tr><th>Predictor</th>${SUMMARY_COLUMNS.map((c) => columnHead(c, c.id)).join("")}</tr>`;
  const body = runs.map((r) => `<tr class="${r.companion_of ? "companion" : ""}">${runCell(r)}${SUMMARY_COLUMNS.map((c) => valueCell(r, c.id, c.fmt, c.tip)).join("")}</tr>`).join("");
  $("#summaryTable").innerHTML = `<table class="grid summary">${head}${body}</table>`;
  const n = (r, id, key = "numerator") => { const m = scopedMetric(r.run_id, S.policy, id); return m && (m.status === "ok" || m.status === "inconsistent") ? fmtInt(m[key]) : '<span class="na">n/a</span>'; };
  $("#summaryCounts").innerHTML = `<table class="grid compact"><tr><th>Predictor</th><th class="num">R reference genes</th><th class="num">Q predicted genes</th><th class="num">loci recovered</th>
    <th class="num">one-to-one exact pairs (TP)</th><th class="num" data-tip="reference genes with a coordinate-exact CDS partner, before one-to-one matching">reference genes with an exact partner</th>
    <th class="num" data-tip="predicted genes with a coordinate-exact CDS partner, before one-to-one matching">predicted genes with an exact partner</th></tr>
    ${runs.map((r) => statusOf(r).ok ? `<tr>${runCell(r)}<td class="num">${n(r, "reference_coding_genes")}</td><td class="num">${n(r, "predicted_coding_genes")}</td><td class="num">${n(r, "locus_recovery_cds_overlap")}</td>
      <td class="num">${n(r, "exact_tp_one_to_one")}</td><td class="num">${n(r, "cds_exact_ref")}</td><td class="num">${n(r, "cds_exact_query")}</td></tr>` : `<tr>${runCell(r)}${naCell(statusOf(r).text, statusOf(r).tip, "num")}<td colspan="5"></td></tr>`).join("")}</table>
    <p class="muted small">Per-gene exact counts can exceed TP when one gene has exact partners in two genes; only one pair per gene counts one-to-one. ${S.scope === WHOLE ? "Values are the comparator's own." : `Scope "${esc(scopeInfo().label)}": reference and predicted genes on other sequences are removed before every value, including the one-to-one matching.`}</p>`;
}
const STRUCTURE_COLUMNS = [
  { id: "cds_intron_chain_best_pair", label: "CDS intron chain identical", den: "multi-segment reference genes" },
  { id: "ref_cds_introns_recovered", label: "Reference CDS introns recovered", den: "unique reference CDS introns" },
  { id: "query_cds_introns_supported", label: "Predicted CDS introns supported", den: null, sub: "% of the run's unique CDS introns" },
  { id: "splits", label: "Reference genes split", den: "reference genes (R)" },
  { id: "merges", label: "Predictions merging genes", den: null, sub: "% of the run's predicted genes (Q)" },
];
function commonDenominator(runs, id) {
  const dens = new Set(runs.filter((r) => statusOf(r).ok).map((r) => scopedMetric(r.run_id, S.policy, id)).filter((m) => m && m.status === "ok").map((m) => m.denominator));
  return dens.size === 1 ? [...dens][0] : null;
}
function renderStructureSummary() {
  const runs = selectedRuns();
  if (!runs.length) { $("#structureSummary").innerHTML = '<p class="na">Select at least one prediction run.</p>'; return; }
  const head = `<tr><th>Predictor</th>${STRUCTURE_COLUMNS.map((c) => {
    const den = c.den ? commonDenominator(runs, c.id) : null;
    return columnHead({ label: c.label, sub: c.den ? `% of ${den === null ? "" : fmtInt(den) + " "}${c.den}` : c.sub }, c.id);
  }).join("")}</tr>`;
  const body = runs.map((r) => `<tr class="${r.companion_of ? "companion" : ""}">${runCell(r)}${STRUCTURE_COLUMNS.map((c) => valueCell(r, c.id, "pct",
    (m, v) => `${fmtInt(m.numerator)} / ${fmtInt(m.denominator)} = ${fmtPct(v, 2)}`).replace('<div class="minibar">', `<div class="cell-sub">${fmtInt((scopedMetric(r.run_id, S.policy, c.id) || {}).numerator)}</div><div class="minibar">`)).join("")}</tr>`).join("");
  $("#structureSummary").innerHTML = `<h2>Structural agreement at the gene and intron level <span class="small muted">(view ${esc(S.policy)}, ${esc(scopeInfo().label)})</span></h2><table class="grid summary">${head}${body}</table>`;
}

// ---------------------------------------------------------------- exploratory splice plots
const LEN_SHORT = ["< 10", "10–30", "30–60", "60–100", "100–200", "200–500", "0.5–1k", "1–2k", "2–5k", "5–10k", "10–20k", "20–50k", "50–100k", "≥ 100k"];
function plotPolicies() { return S.meta.policies.filter((p) => p.reference_transcript_selection !== "all"); }
function seg(el, options, value, onPick) {
  el.innerHTML = options.map(([v, label, tip]) => `<button data-v="${esc(v)}" class="${v === value ? "active" : ""}" ${tip ? `data-tip="${esc(tip)}"` : ""}>${esc(label)}</button>`).join("");
  $$("button", el).forEach((b) => (b.onclick = () => onPick(b.dataset.v)));
}
async function renderPlots() {
  const pols = plotPolicies();
  if (!pols.length) { $("#plotStatus").innerHTML = '<p class="na">Not available: every policy evaluates all reference isoforms; these plots need one reference transcript per gene.</p>'; return; }
  if (!S.plotPolicy || !pols.some((p) => p.id === S.plotPolicy)) S.plotPolicy = (pols.find((p) => p.role !== "supplementary") || pols[0]).id;
  seg($("#plotPolicy"), pols.map((p) => [p.id, `View ${p.id}`, `${p.label}${p.role === "supplementary" ? " (supplementary)" : ""}`]), S.plotPolicy, (v) => { S.plotPolicy = v; renderPlots(); });
  seg($("#posCoord"), [["cds_span", "genomic CDS span", "position of the exon–intron junction along the genomic span from coding start to coding end (introns included)"],
    ["coding", "coding sequence", "fraction of coding bases upstream of the junction"]], S.posCoord, (v) => { S.posCoord = v; drawPlots(); });
  seg($("#posMin"), [["0", "all CDS introns", "every reference CDS intron"],
    ["60", "introns ≥ 60 bp", "leave out reference CDS introns shorter than 60 bp: mostly gaps of a few bases that predictors rarely model (see the intron-length plot)"]],
    String(S.posMin), (v) => { S.posMin = Number(v); drawPlots(); });
  seg($("#lenMeasure"), [["exact", "complete intron", "both boundaries match: an identical CDS intron"], ["donor", "donor site", "the 5' splice site matches"],
    ["acceptor", "acceptor site", "the 3' splice site matches"]], S.lenMeasure, (v) => { S.lenMeasure = v; drawPlots(); });
  const p = policyInfo(S.plotPolicy);
  const rec = (S.ov.analyses || []).find((a) => a.policy_id === S.plotPolicy && a.record?.support_rule)?.record;
  $("#plotNote").innerHTML = `Plots use <b>view ${esc(p.id)}</b> (${esc(p.label)}) with one reference transcript per gene, so each reference intron counts once. The headline policy selector does not change them; view A is not offered because a gene with several isoforms would contribute shared introns several times.
    Support rule: ${esc(rec?.support_rule || "exact same-sequence, same-strand CDS intron or splice-site coordinates in any evaluated predicted transcript")}. Sites in reference genes whose locus was not recovered are <b>missing-locus sites</b>: counted separately and excluded from the rates, so the rates describe structural errors at detected loci. Scope: ${esc(scopeInfo().label)}.`;
  const key = `${S.ref}|${S.plotPolicy}|${S.scope}`;
  if (!S.splice[key]) {
    try { S.splice[key] = await api("/api/splice", { reference: S.ref, policy: S.plotPolicy, scope: S.scope }); }
    catch (e) { S.splice[key] = { error: e.message.replace(/^QueryError: /, "") }; }
  }
  drawPlots();
}
function analysisIssues(runs) {
  const out = [];
  for (const r of runs) {
    const a = (S.ov.analyses || []).find((x) => x.run_id === r.run_id && x.policy_id === S.plotPolicy);
    if (!a) out.push(`${runDisplay(r)}: not analysed`);
    else if (a.status !== "ok") out.push(`${runDisplay(r)}: ${a.status}${a.reason ? " — " + a.reason : ""}`);
  }
  return out;
}
function drawPlots() {
  const data = S.splice[`${S.ref}|${S.plotPolicy}|${S.scope}`];
  if (!data) return;
  const runs = selectedRuns();
  for (const id of ["#positionChart", "#positionTable", "#lengthChart", "#lengthTable", "#positionNote", "#lengthNote", "#spliceDownloads"]) $(id).innerHTML = "";
  if (data.error) { $("#plotStatus").innerHTML = `<p class="na">Not available: ${esc(data.error)}</p>`; return; }
  const issues = analysisIssues(runs);
  const ok = runs.filter((r) => data.sites.some((x) => x.run_id === r.run_id));
  $("#plotStatus").innerHTML = issues.length ? `<div class="issue warning">${issues.map(esc).join("<br>")}</div>` : "";
  if (!ok.length) { $("#positionChart").innerHTML = '<p class="na">No analysed results for the selected runs.</p>'; return; }
  drawPositionChart(data, ok); drawLengthChart(data, ok); writeHash();
  $("#spliceDownloads").innerHTML = `<span class="small muted">Binned counts (TSV, view ${esc(S.plotPolicy)}, ${esc(scopeInfo().label)}):</span>
    <a class="btn" href="${esc(apiUrl("/api/splice", { reference: S.ref, policy: S.plotPolicy, scope: S.scope, format: "csv", table: "sites" }))}">splice sites by position</a>
    <a class="btn" href="${esc(apiUrl("/api/splice", { reference: S.ref, policy: S.plotPolicy, scope: S.scope, format: "csv", table: "lengths" }))}">introns by length</a>`;
}
function niceMax(v) { return Math.min(1, Math.max(0.1, Math.ceil(v * 10 + 1e-9) / 10)); }
function lineChart({ width, height, labels, series, yMax, yTitle, xTitle, footer, edgeLabels }) {
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
  for (const s of series) {
    let path = "";
    s.points.forEach((p, i) => { if (p.v === null) { path += " "; return; } path += `${path.endsWith(" ") || !path ? "M" : "L"}${x(i).toFixed(1)},${y(p.v).toFixed(1)}`; });
    out += `<path d="${path.replace(/ /g, "")}" fill="none" stroke="${s.color}" stroke-width="2" ${s.dash ? 'stroke-dasharray="5 3"' : ""}/>`;
    s.points.forEach((p, i) => { if (p.v !== null) out += `<circle cx="${x(i)}" cy="${y(p.v)}" r="3.5" fill="${s.color}" data-tip="${esc(p.tip)}"/>`; });
  }
  return svg(width, height, out);
}
function legendHtml(runs) {
  return `<div class="legend">${runs.map((r) => `<span><span class="swatch" style="background:${r.color}"></span> ${esc(runDisplay(r))}${r.companion_of ? " (supplementary)" : ""}</span>`).join("")}</div>`;
}
function rate(num, den) { return den ? num / den : null; }
function drawPositionChart(data, runs) {
  const el = $("#positionChart");
  const W = Math.min(1100, el.clientWidth || 1000), half = W >= 760 ? Math.floor((W - 16) / 2) : W;
  const labels = Array.from({ length: 10 }, (_, i) => `${(i / 10).toFixed(1)}–${((i + 1) / 10).toFixed(1)}`);
  const rows = (run, kind) => data.sites.filter((x) => x.run_id === run && (x.min_intron_length ?? 0) === S.posMin && x.coordinate === S.posCoord && x.kind === kind).sort((a, b) => a.bin - b.bin);
  let max = 0;
  for (const r of runs) for (const k of ["donor", "acceptor"]) for (const b of rows(r.run_id, k)) { const v = rate(b.assessed - b.supported, b.assessed); if (v !== null) max = Math.max(max, v); }
  const yMax = niceMax(max);
  const panel = (kind) => lineChart({
    width: half, height: 250, labels, yMax, yTitle: "unsupported (% of assessed)", xTitle: `${kind} sites · position along the reference CDS (0 = coding start, 1 = coding end)`,
    edgeLabels: (i) => (i / 10).toFixed(1),
    series: runs.map((r) => ({ color: r.color, dash: !!r.companion_of, points: rows(r.run_id, kind).map((b) => ({
      v: rate(b.assessed - b.supported, b.assessed),
      tip: `${runDisplay(r)} · ${kind} sites, position ${labels[b.bin]}\nunsupported ${fmtInt(b.assessed - b.supported)} of ${fmtInt(b.assessed)} assessed (${b.assessed ? fmtPct((b.assessed - b.supported) / b.assessed) : "not available"})\n${fmtInt(b.missing_locus)} further sites in genes whose locus was not recovered (excluded)` })) })),
  });
  el.innerHTML = legendHtml(runs) + `<div class="panels">${panel("donor")}${panel("acceptor")}</div>`;
  const miss = runs.map((r) => {
    const b = rows(r.run_id, "donor"); const a = b.reduce((t, x) => t + x.assessed, 0), m = b.reduce((t, x) => t + x.missing_locus, 0);
    return `${esc(runDisplay(r))} ${fmtInt(m)} of ${fmtInt(a + m)} (${a + m ? fmtPct(m / (a + m)) : "—"})`;
  });
  $("#positionNote").innerHTML = `A reference donor (acceptor) is <b>unsupported</b> when no evaluated predicted transcript on the same sequence and strand has a CDS intron with that donor (acceptor) coordinate. Genes are oriented in the direction of transcription, minus-strand genes included.
    Each bin's rate is over the sites assessed in it (counts in the table); empty bins are not drawn. This is a distribution by position, not a cumulative count.
    ${S.posMin ? `Reference CDS introns shorter than ${S.posMin} bp are left out (toggle above).` : "All reference CDS introns are included, also gaps of a few bases that predictors rarely model; use “introns ≥ 60 bp” to leave them out."}
    <br><b>Missing-locus sites</b> (excluded above; donor sites, same for acceptors): ${miss.join(" · ")}.`;
  const head = `<tr><th>Predictor</th><th>site</th>${labels.map((l) => `<th class="num">${l}</th>`).join("")}<th class="num">missing-locus sites</th></tr>`;
  const body = runs.flatMap((r) => ["donor", "acceptor"].map((k) => { const b = rows(r.run_id, k);
    return `<tr><td>${esc(runDisplay(r))}</td><td>${k}</td>${b.map((x) => x.assessed ? `<td class="num" data-tip="${esc(`${fmtInt(x.supported)} supported, ${fmtInt(x.assessed - x.supported)} unsupported`)}">${fmtPct((x.assessed - x.supported) / x.assessed)}<div class="cell-sub">${fmtInt(x.assessed - x.supported)} / ${fmtInt(x.assessed)}</div></td>` : '<td class="num"><span class="na">not available</span><div class="cell-sub">0 assessed</div></td>').join("")}
      <td class="num">${fmtInt(b.reduce((t, x) => t + x.missing_locus, 0))}</td></tr>`; })).join("");
  $("#positionTable").innerHTML = `<table class="grid compact">${head}${body}</table><p class="small muted">Cells: unsupported / assessed sites at recovered loci (rate). Coordinate: ${S.posCoord === "coding" ? "coding sequence" : "genomic CDS span"}; ${S.posMin ? `reference CDS introns ≥ ${S.posMin} bp` : "all reference CDS introns"}.</p>`;
}
function drawLengthChart(data, runs) {
  const el = $("#lengthChart");
  const W = Math.min(1100, el.clientWidth || 1000);
  const rows = (run) => data.lengths.filter((x) => x.run_id === run).sort((a, b) => a.bin - b.bin);
  const first = rows(runs[0].run_id);
  const refTotals = first.map((b) => b.eligible + b.missing_locus);
  const consistent = runs.every((r) => rows(r.run_id).every((b, i) => b.eligible + b.missing_locus === refTotals[i]));
  const m = S.lenMeasure, word = { exact: "complete introns recovered", donor: "donor sites supported", acceptor: "acceptor sites supported" }[m];
  el.innerHTML = legendHtml(runs) + lineChart({
    width: W, height: 290, labels: first.map((b) => LEN_SHORT[b.bin] || b.label), yMax: 1, yTitle: `${word} (% of eligible)`,
    xTitle: "reference CDS intron length (bp, fixed bins, same for every run)", footer: refTotals.map((n) => `n=${fmtInt(n)}`),
    series: runs.map((r) => ({ color: r.color, dash: !!r.companion_of, points: rows(r.run_id).map((b) => ({
      v: rate(b[m], b.eligible),
      tip: `${runDisplay(r)} · introns ${b.label}\n${word}: ${fmtInt(b[m])} of ${fmtInt(b.eligible)} at recovered loci (${b.eligible ? fmtPct(b[m] / b.eligible) : "not available"})\ncomplete ${fmtInt(b.exact)} · donor ${fmtInt(b.donor)} · acceptor ${fmtInt(b.acceptor)}\n${fmtInt(b.missing_locus)} introns in genes whose locus was not recovered (excluded)` })) })),
  });
  $("#lengthNote").innerHTML = `<b>Complete intron</b>: an evaluated predicted transcript has the identical CDS intron (both boundaries). <b>Donor / acceptor</b>: only that splice site matches. n = reference CDS introns in each bin (one transcript per gene${consistent ? "" : "; differs between runs — check the analysis status"}); rates are over introns of genes whose locus was recovered.
    These are associations, not causes: short- and long-intron bins hold different genes and different numbers of introns.`;
  const head = `<tr><th>Predictor</th>${first.map((b) => `<th class="num">${esc(b.label)}</th>`).join("")}</tr>`;
  const body = runs.map((r) => `<tr><td>${esc(runDisplay(r))}</td>${rows(r.run_id).map((b) => b.eligible
    ? `<td class="num" data-tip="${esc(`complete ${fmtInt(b.exact)}, donor ${fmtInt(b.donor)}, acceptor ${fmtInt(b.acceptor)} of ${fmtInt(b.eligible)}; ${fmtInt(b.missing_locus)} in missing loci`)}">${fmtPct(b[m] / b.eligible)}<div class="cell-sub">${fmtInt(b[m])} / ${fmtInt(b.eligible)}</div><div class="cell-sub">missing loci ${fmtInt(b.missing_locus)}</div></td>`
    : `<td class="num"><span class="na">not available</span><div class="cell-sub">0 eligible · missing loci ${fmtInt(b.missing_locus)}</div></td>`).join("")}</tr>`).join("");
  $("#lengthTable").innerHTML = `<table class="grid compact">${head}${body}</table><p class="small muted">Cells: ${esc(word)} / eligible introns at recovered loci; missing-locus introns are excluded from the rate.</p>`;
}

// ---------------------------------------------------------------- charts (SVG)
function svg(w, h, inner) { return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img">${inner}</svg>`; }
const CHAR_W = 6.4; // approximate width of a 12px label character
function chartLabel(r, maxChars) {
  const name = runDisplay(r);
  return name.length > maxChars ? `${esc(name.slice(0, maxChars - 1))}…` : esc(name);
}
function labelColumn(runs, W, indent) {
  const longest = Math.max(...runs.map((r) => runDisplay(r).length));
  return Math.round(Math.min(Math.max(210, indent + longest * CHAR_W + 16), W * 0.45));
}
function renderPolicyChart() {
  const id = $("#policyMetric").value; const runs = selectedRuns(); const pols = S.meta.policies;
  const el = $("#policyChart");
  if (!runs.length) { el.innerHTML = ""; return; }
  const W = Math.min(980, el.clientWidth || 980), L = labelColumn(runs, W, 14), bh = 16, gap = 4, groupGap = 14;
  let y = 8, out = "";
  const x = (v) => L + v * (W - L - 130);
  for (const p of pols) {
    out += `<text x="4" y="${y + 12}" font-weight="600">${esc(p.id)} · ${esc(p.label)}</text>`;
    y += 18;
    for (const r of runs) {
      const res = result(r.run_id, p.id); const m = metric(r.run_id, p.id, id);
      out += `<text x="14" y="${y + 12}" class="lbl-muted">${chartLabel(r, Math.floor((L - 30) / CHAR_W))}</text>`;
      if (!res || res.status !== "complete") {
        out += `<text x="${L}" y="${y + 12}" fill="var(--na)" font-style="italic">${esc(res ? RESULT_LABEL[res.status] : "no result")}${res?.reason ? " — " + esc(res.reason.slice(0, 80)) : ""}</text>`;
      } else if (!m || m.status !== "ok" || !m.denominator) {
        out += `<text x="${L}" y="${y + 12}" fill="var(--na)" font-style="italic">${esc(m ? METRIC_STATUS[m.status] || m.status : "n/a")} — ${esc((m?.note || "").slice(0, 90))}</text>`;
      } else {
        const v = m.value ?? m.numerator / m.denominator;
        out += `<rect x="${L}" y="${y + 2}" width="${Math.max(1, x(v) - L)}" height="${bh - 4}" rx="3" fill="${r.color}" data-tip="${esc(`${runLabel(r, false)}, policy ${p.id}\n${fmtInt(m.numerator)} / ${fmtInt(m.denominator)} = ${fmtPct(v, 2)}`)}"></rect>`;
        const shown = S.dict[id]?.unit === "rate" ? v.toFixed(4) : fmtPct(v);
        out += `<text x="${x(v) + 6}" y="${y + 12}">${shown} <tspan fill="var(--muted)">(${fmtInt(m.numerator)} / ${fmtInt(m.denominator)})</tspan></text>`;
      }
      y += bh + gap;
    }
    y += groupGap;
  }
  el.innerHTML = svg(W, y, out);
}
function stacked(el, title, parts, totalId) {
  const runs = selectedRuns().filter((r) => result(r.run_id, S.policy)?.status === "complete");
  if (!runs.length) { el.innerHTML = ""; return; }
  const W = Math.min(980, el.clientWidth || 980), L = labelColumn(runs, W, 4), bh = 18;
  let out = `<text x="4" y="14" font-weight="600">${esc(title)}</text>`;
  // legend: flow layout on its own lines below the title
  let lx = 4, ly = 24;
  for (const p of parts) {
    const w = 22 + p.label.length * 6.4;
    if (lx + w > W - 4) { lx = 4; ly += 16; }
    out += `<rect x="${lx}" y="${ly}" width="10" height="10" fill="${p.color}" rx="2"></rect><text x="${lx + 14}" y="${ly + 9}">${esc(p.label)}</text>`;
    lx += w + 14;
  }
  let y = ly + 22;
  for (const r of runs) {
    const total = metric(r.run_id, S.policy, totalId)?.numerator;
    out += `<text x="4" y="${y + 13}">${chartLabel(r, Math.floor((L - 20) / CHAR_W))}</text>`;
    let x = L;
    for (const p of parts) {
      const n = p.value(r);
      if (n === null || n === undefined || !total) continue;
      const w = (n / total) * (W - L - 20);
      out += `<rect x="${x}" y="${y}" width="${Math.max(0, w - 2)}" height="${bh}" fill="${p.color}" rx="2" data-tip="${esc(`${runLabel(r, false)}\n${p.label}: ${fmtInt(n)} / ${fmtInt(total)} (${fmtPct(n / total)})`)}"></rect>`;
      x += w;
    }
    y += bh + 6;
  }
  el.innerHTML = svg(W, y + 4, out);
}
function renderPartialChart() {
  const v = (id) => (r) => { const m = metric(r.run_id, S.policy, id); return m && m.status === "ok" ? m.numerator : null; };
  stacked($("#partialChart"), "Reference genes by CDS outcome (sums to R)", [
    { label: "coordinate-exact", color: "#1e6b3f", value: v("cds_exact_ref") },
    { label: "ends differ", color: "#c98500", value: v("partial_same_chain_terminal_diff") },
    { label: "splice differs", color: "#eda100", value: v("partial_ge08_diff_chain") },
    { label: "< 0.8 overlap", color: "#f3cd75", value: v("partial_lt08") },
    { label: "no CDS overlap", color: "#c0453f", value: v("ref_genes_without_cds_overlap") },
  ], "reference_coding_genes");
}
function renderNovelChart() {
  stacked($("#novelChart"), "Unmatched predictions by reference context (share of all predictions Q)", NOVEL.map((id, i) => ({
    label: (S.dict[id]?.label || id).replace("Novel: ", ""), color: NOVEL_COLORS[i],
    value: (r) => { const m = metric(r.run_id, S.policy, id); return m && m.status === "ok" ? m.numerator : null; } })), "predicted_coding_genes");
}

// ---------------------------------------------------------------- gene table
function setupGeneFilters() {
  const runs = refRuns().filter((r) => S.meta.runs);
  $("#focusRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(runDisplay(r))}</option>`).join("");
  $("#compareRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(runDisplay(r))}</option>`).join("");
  if (runs[1]) $("#compareRun").value = runs[1].run_id;
  $("#statusChecks").innerHTML = Object.entries(STATUS).map(([k, [lab, t]]) => `<label data-tip="${esc(t)}"><input type="checkbox" value="${k}"> <span class="status st-${k}">${esc(lab)}</span></label>`).join("");
  api("/api/chromosomes", { reference: S.ref }).then((chroms) => {
    $("#chrom").innerHTML = '<option value="">all sequences</option>' + chroms.map((c) => `<option>${esc(c)}</option>`).join("");
  });
  $("#geneFilters").onsubmit = (e) => { e.preventDefault(); S.genePage = 0; loadGenes(); };
  $("#resetFilters").onclick = () => { $("#geneFilters").reset(); S.genePage = 0; loadGenes(); };
  $("#dlGenes").onclick = () => (location.href = apiUrl("/api/genes", { ...geneParams(), format: "csv", runs: selectedRuns().map((r) => r.run_id).join(",") }));
}
function geneParams() {
  const statuses = $$("#statusChecks input:checked").map((c) => c.value).join(",");
  return { reference: S.ref, policy: S.policy, q: $("#q").value, chrom: $("#chrom").value, multi: $("#multi").value,
    focus_run: $("#focusRun").value, status: statuses, split: $("#splitOnly").checked ? "1" : "", contrast: $("#contrast").value,
    compare_run: $("#compareRun").value, limit: S.pageSize, offset: S.genePage * S.pageSize };
}
async function loadGenes() {
  S.genesLoaded = true;
  const token = (S.geneToken = (S.geneToken || 0) + 1);
  try {
    const data = await api("/api/genes", geneParams());
    if (token !== S.geneToken) return; // a newer request superseded this one
    S.lastGenes = data; renderGeneTable(data);
  } catch (e) { showError($("#geneTable"), e); }
}
function statusChip(o) {
  if (!o) return '<span class="status st-none" data-tip="no outcome for this run (failed, not run, unavailable or gene not evaluated)">—</span>';
  const [lab, t] = STATUS[o.cds_status] || [o.cds_status, ""];
  const extra = [`comparator: ${o.classification} / CDS ${o.classification_cds}`, `CDS overlap ${Number(o.cds_overlap).toFixed(3)}`,
    o.cds_matched_id ? `partner ${o.cds_matched_id}` : "", o.one_to_one_partner ? `one-to-one exact pair: ${o.one_to_one_partner}` : (o.cds_exact ? "exact, but not in the one-to-one matching (partner used by another gene)" : ""),
    o.counterpart_count >= 2 ? `split: ${o.counterpart_ids}` : "",
    o.any_pair_chain === null || o.any_pair_chain === undefined ? "" : `any-pair intron chain (independent): ${o.any_pair_chain ? "yes" : "no"}`].filter(Boolean).join("\n");
  return `<span class="status st-${esc(o.cds_status)}" data-tip="${esc(t + "\n" + extra)}">${esc(lab)}${o.counterpart_count >= 2 ? " ⇉" : ""}${o.cds_exact && !o.one_to_one_partner && o.one_to_one_partner !== undefined ? " ≠1:1" : ""}</span>`;
}
function renderGeneTable(data) {
  if (!data) return;
  const runs = selectedRuns();
  $("#geneCount").textContent = `${fmtInt(data.total)} evaluated reference genes match (policy ${S.policy}); showing ${data.offset + 1}–${data.offset + data.rows.length}. ⇉ = split.`;
  const head = `<tr><th>Reference gene</th><th>Name (label)</th><th>Location</th><th>CDS</th>${runs.map((r) => `<th>${runLabel(r)}</th>`).join("")}</tr>`;
  const body = data.rows.map((g) => `<tr class="click ${S.selectedGene === g.gene_id ? "selected" : ""}" data-gene="${esc(g.gene_id)}">
      <td class="mono">${esc(g.gene_id)}${g.original_gene_id && g.original_gene_id !== g.gene_id ? `<div class="cell-sub">file ID ${esc(g.original_gene_id)}</div>` : ""}</td>
      <td>${esc(g.name || "")}</td><td class="mono">${esc(g.chrom)}:${fmtInt(g.start)}-${fmtInt(g.end)} ${esc(g.strand)}</td>
      <td class="small">${g.multi_segment ? "multi" : "single"}-segment · ${g.evaluated_transcript_count} isoform${g.evaluated_transcript_count === 1 ? "" : "s"}${g.canonical_fallback ? ' <span class="tag warn" data-tip="no Ensembl_canonical tag: longest CDS used">fallback</span>' : ""}</td>
      ${runs.map((r) => `<td>${statusChip(g.outcomes[r.run_id])}</td>`).join("")}</tr>`).join("");
  $("#geneTable").innerHTML = `<table class="grid">${head}${body}</table>`;
  $$("#geneTable tr.click").forEach((tr) => (tr.onclick = () => {
    const g = data.rows.find((x) => x.gene_id === tr.dataset.gene); S.selectedGene = g.gene_id;
    $$("#geneTable tr").forEach((x) => x.classList.toggle("selected", x === tr)); openLocus(g.chrom, g.start, g.end, g.gene_id);
  }));
  pager($("#genePager"), data.total, S.genePage, (p) => { S.genePage = p; loadGenes(); });
}
function pager(el, total, page, go) {
  const pages = Math.max(1, Math.ceil(total / S.pageSize));
  el.innerHTML = `<button class="btn" ${page <= 0 ? "disabled" : ""} data-p="${page - 1}">‹ previous</button><span class="small">page ${page + 1} of ${pages}</span>
    <button class="btn" ${page >= pages - 1 ? "disabled" : ""} data-p="${page + 1}">next ›</button>`;
  $$("button", el).forEach((b) => (b.onclick = () => go(Number(b.dataset.p))));
}

// ---------------------------------------------------------------- predictions table
function setupPredFilters() {
  const runs = refRuns();
  $("#predRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(runDisplay(r))}</option>`).join("");
  $("#pnovel").innerHTML = '<option value="">any Novel context</option>' + NOVEL.map((id) => `<option value="${id.replace("novel_", "")}">${esc(S.dict[id]?.label || id)}</option>`).join("");
  $("#predFilters").onsubmit = (e) => { e.preventDefault(); S.predPage = 0; loadPredictions(); };
  $("#dlPreds").onclick = () => (location.href = apiUrl("/api/predictions", { ...predParams(), format: "csv" }));
}
function predParams() {
  return { reference: S.ref, policy: S.policy, run: $("#predRun").value, q: $("#pq").value, classification: $("#pclass").value,
    novel_category: $("#pnovel").value, exact: $("#pexact").value, merge: $("#pmerge").checked ? "1" : "", limit: S.pageSize, offset: S.predPage * S.pageSize };
}
async function loadPredictions() {
  S.predsLoaded = true;
  const res = result($("#predRun").value, S.policy);
  if (!res || res.status !== "complete") { $("#predTable").innerHTML = `<p class="na">No complete result for this run and policy (${esc(res ? RESULT_LABEL[res.status] : "no result")}).</p>`; $("#predCount").textContent = ""; $("#predPager").innerHTML = ""; return; }
  const token = (S.predToken = (S.predToken || 0) + 1);
  try {
    const params = predParams();
    const data = await api("/api/predictions", params);
    if (token !== S.predToken) return; // superseded
    const run = runById(params.run);
    $("#predCount").textContent = `${fmtInt(data.total)} predictions match; showing ${data.offset + 1}–${data.offset + data.rows.length}.`;
    const head = "<tr><th>Prediction (comparison ID)</th><th>Location</th><th>Classification</th><th>CDS exact</th><th>Reference partner</th><th>Novel context</th><th>Counterparts</th></tr>";
    const body = data.rows.map((q) => `<tr class="click" data-chrom="${esc(q.chrom)}" data-s="${q.start}" data-e="${q.end}">
      <td class="mono">${esc(q.gene_id)}${q.original_gene_id !== q.gene_id ? `<div class="cell-sub">file ID ${esc(q.original_gene_id)} (namespaced: ID reused on several sequences)</div>` : ""}</td>
      <td class="mono">${esc(q.chrom)}:${fmtInt(q.start)}-${fmtInt(q.end)} ${esc(q.strand)}</td>
      <td>${esc(q.classification)}<div class="cell-sub">CDS ${esc(q.classification_cds)}</div></td>
      <td>${q.cds_exact ? '<span class="status st-exact">exact</span>' : '<span class="small muted">no</span>'}</td>
      <td class="mono">${esc(q.cds_matched_id || "")}<div class="cell-sub">${esc(q.best_cds_ref_tx || "")}</div></td>
      <td class="small">${esc((q.novel_category || "").replaceAll("_", " "))}</td>
      <td class="small">${q.counterpart_count >= 2 ? `<span class="tag warn">merge</span> ${esc(q.counterpart_ids)}` : esc(q.counterpart_ids || "")}</td></tr>`).join("");
    $("#predTable").innerHTML = `<table class="grid">${head}${body}</table>`;
    $$("#predTable tr.click").forEach((tr) => (tr.onclick = () => { selectTab("genes"); S.runs.add(run.run_id); openLocus(tr.dataset.chrom, Number(tr.dataset.s), Number(tr.dataset.e), null); }));
    pager($("#predPager"), data.total, S.predPage, (p) => { S.predPage = p; loadPredictions(); });
  } catch (e) { showError($("#predTable"), e); }
}

// ---------------------------------------------------------------- locus browser
function openLocus(chrom, start, end, geneId) {
  const ref = S.meta.references.find((r) => r.id === S.ref);
  $("#locusPanel").hidden = false;
  if (!ref.models_available) {
    $("#locusSvg").innerHTML = '<p class="na">Locus view needs the prepared reference annotation and the prediction files; this dataset was built without them (summary-only). Gene outcomes above come from comparison_details.tsv.</p>';
    return;
  }
  const pad = Math.max(500, Math.round((end - start) * 0.15));
  S.locus = { chrom, start: Math.max(1, start - pad), end: end + pad, gene: geneId };
  drawLocus();
  $("#locusPanel").scrollIntoView({ behavior: "smooth", block: "start" });
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
function pack(items, key = (t) => [t.start, t.end], padBp = 0) {
  const rows = [];
  for (const it of items) {
    const [s, e] = key(it);
    let row = rows.findIndex((last) => last < s - padBp);
    if (row < 0) { rows.push(e); row = rows.length - 1; } else rows[row] = e;
    it._row = row;
  }
  return rows.length;
}
async function drawLocus() {
  const L = S.locus; const runs = selectedRuns();
  $("#locusInput").value = `${L.chrom}:${L.start}-${L.end}`;
  let data;
  try {
    data = await api("/api/locus", { reference: S.ref, chrom: L.chrom, start: L.start, end: L.end, policy: S.policy,
      runs: runs.map((r) => r.run_id).join(","), gene_id: L.gene || "", isoforms: $("#isoforms").value });
  } catch (e) { showError($("#locusSvg"), e); return; }
  L.start = data.start; L.end = data.end;
  $("#locusTitle").textContent = `Locus ${data.chrom}:${fmtInt(data.start)}-${fmtInt(data.end)} (${fmtInt(data.end - data.start + 1)} bp)`;
  const notes = [];
  if (data.window_clamped) notes.push(`window limited to ${fmtInt(data.limits.max_window)} bp`);
  if (data.reference_transcripts_truncated) notes.push(`reference transcripts limited to ${data.limits.track_limit}`);
  if (data.reference_transcripts_hidden) notes.push(`${data.reference_transcripts_hidden} other reference isoforms hidden (change the isoform selector)`);
  if (data.reference_genes_truncated) notes.push(`reference gene context limited to ${data.limits.context_limit}`);
  for (const [run, t] of Object.entries(data.runs)) if (t.truncated) notes.push(`${run}: limited to ${data.limits.track_limit} transcripts`);
  $("#locusNotes").textContent = notes.join(" · ");
  $("#locusLegend").innerHTML = `<span><svg width="34" height="14"><line x1="0" y1="7" x2="34" y2="7" stroke="var(--ref)"/><rect x="4" y="4" width="10" height="6" fill="none" stroke="var(--ref)"/><rect x="18" y="1" width="12" height="12" fill="var(--ref)"/></svg> UTR exon (thin, outlined) / CDS (thick, filled) / intron (line)</span>
    <span>★ Ensembl_canonical tag · ◆ compared isoform (best CDS pair) · faded = not evaluated under this policy</span>
    <span><svg width="22" height="10"><rect width="22" height="10" fill="url(#hatch)" stroke="var(--na)"/></svg> reference gene not evaluated (e.g. pseudogene, non-coding)</span>`;

  const W = Math.max(760, $("#locusSvg").clientWidth || 1100), LBL = 230, R = W - 12, rowH = 20;
  const x = (p) => LBL + ((p - data.start) / (data.end - data.start + 1)) * (R - LBL);
  const clip = (p) => Math.min(R, Math.max(LBL, x(p)));
  let y = 6, out = `<defs><pattern id="hatch" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="6" stroke="var(--na)" stroke-width="2"/></pattern></defs>`;
  // axis
  const span = data.end - data.start + 1, step = Math.pow(10, Math.floor(Math.log10(span / 5)));
  const tickStep = span / step > 10 ? step * 2 : step;
  out += `<g class="axis"><line x1="${LBL}" y1="${y + 16}" x2="${R}" y2="${y + 16}"/>`;
  for (let t = Math.ceil(data.start / tickStep) * tickStep; t <= data.end; t += tickStep)
    out += `<line x1="${x(t)}" y1="${y + 12}" x2="${x(t)}" y2="${y + 20}"/><text x="${x(t)}" y="${y + 9}" text-anchor="middle" class="lbl-muted">${fmtInt(t)}</text>`;
  out += `</g>`; y += 30;
  // focus band
  const focusGene = data.reference_genes.find((g) => g.gene_id === L.gene);
  const bandTop = y;
  // gene context
  const ctx = data.reference_genes;
  const nrow = pack(ctx, (g) => [x(g.start), x(g.start) + Math.max(x(g.end) - x(g.start), 8 * (g.name || g.gene_id).length)], 4);
  out += `<text x="4" y="${y + 12}" font-weight="600">Reference genes</text>`;
  for (const g of ctx) {
    const gy = y + g._row * 16, ev = g.evaluated;
    const fill = ev ? "var(--ref)" : "url(#hatch)";
    out += `<g data-tip="${esc(`${g.gene_id}${g.name ? " (" + g.name + ")" : ""}\n${g.biotype} ${g.chrom}:${g.start}-${g.end} ${g.strand}\n${ev ? "evaluated" : "not evaluated (filtered by biotype/scope)"}`)}">
      <rect x="${clip(g.start)}" y="${gy + 4}" width="${Math.max(2, clip(g.end) - clip(g.start))}" height="6" fill="${fill}" stroke="${ev ? "none" : "var(--na)"}" opacity="${ev ? 0.85 : 1}"/>
      <text x="${clip(g.start)}" y="${gy + 3}" class="${ev ? "" : "lbl-muted"}" font-size="10">${esc(g.name || g.gene_id)} ${g.strand === "-" ? "◀" : "▶"}</text></g>`;
  }
  y += Math.max(1, nrow) * 16 + 12;
  // reference transcripts
  out += `<text x="4" y="${y + 12}" font-weight="600">Reference transcripts (${esc(S.policy)})</text>`; y += 18;
  for (const t of data.reference_transcripts) {
    out += transcriptGlyph(t, y, x, clip, "var(--ref)", t.evaluated ? 1 : 0.35,
      `${t.canonical ? "★ " : ""}${t.partner_of_runs.length ? "◆ " : ""}${t.transcript_id}`,
      `${t.transcript_id} (${t.biotype || "?"})\ngene ${t.gene_id}\n${t.evaluated ? "evaluated" : "not evaluated"} under policy ${S.policy}${t.tags ? "\ntags: " + t.tags : ""}${t.partner_of_runs.length ? "\nbest CDS pair for: " + t.partner_of_runs.join(", ") : ""}`,
      t.partner_of_runs.map((r) => runById(r)?.color));
    y += rowH;
  }
  if (!data.reference_transcripts.length) { out += `<text x="${LBL}" y="${y + 12}" class="lbl-muted">no reference transcripts in this window for the selected isoform view</text>`; y += rowH; }
  // runs
  for (const r of runs) {
    const tr = data.runs[r.run_id];
    y += 6;
    out += `<rect x="4" y="${y + 3}" width="10" height="10" fill="${r.color}" rx="2"/><text x="18" y="${y + 12}" font-weight="600">${esc(runDisplay(r))} <tspan class="lbl-muted" font-weight="400">${harmonised(r) ? "evaluated (stop-harmonised) coordinates" : ""}</tspan></text>`;
    y += 18;
    const modelError = S.meta.model_errors?.[r.run_id];
    if (modelError || !tr || !tr.transcripts.length) {
      const msg = modelError ? `prediction models unavailable — ${modelError.slice(0, 110)}` : tr ? "no predictions in this window" : "prediction models unavailable";
      out += `<text x="${LBL}" y="${y + 12}" class="lbl-muted">${esc(msg)}</text>`; y += rowH; continue;
    }
    for (const t of tr.transcripts) {
      const o = t.outcome, st = t.stop_codon;
      const stopLabel = st ? (st.status === "extended" ? " · stop +3" : ` · stop ${st.status}`) : "";
      const lab = `${t.gene_id}${o && o.cds_exact ? " · CDS exact" : ""}${stopLabel}${o ? " · " + o.classification : ""}`;
      const stopTip = st ? `\nstop codon (--add-stop-codon): ${st.status} — ${st.reason}${st.last_codon ? "\nlast codon " + st.last_codon : ""}${st.next_codon ? ", next codon " + st.next_codon : ""}\nterminal CDS ${st.terminal_cds_start}-${st.terminal_cds_end}${st.status === "extended" ? " → " + st.new_terminal_cds_start + "-" + st.new_terminal_cds_end : " (unchanged)"}` : "";
      out += transcriptGlyph(t, y, x, clip, r.color, 1, lab, `${t.transcript_id}\ngene ${t.gene_id}\n${o ? `comparator: ${o.classification} / CDS ${o.classification_cds}${o.novel_category ? "\nNovel: " + o.novel_category : ""}${o.counterpart_count >= 2 ? "\nmerge (" + o.counterpart_count + " reference counterparts)" : ""}` : "no outcome for this policy"}${stopTip}`, []);
      y += rowH;
    }
  }
  if (focusGene) out = `<rect x="${clip(focusGene.start)}" y="${bandTop}" width="${Math.max(2, clip(focusGene.end) - clip(focusGene.start))}" height="${y - bandTop}" fill="var(--info-bg)" opacity="0.6"/>` + out;
  $("#locusSvg").innerHTML = svg(W, y + 8, out);
  if (L.gene) loadGeneDetail(L.gene); else $("#geneDetail").innerHTML = "";
}
function transcriptGlyph(t, y, x, clip, color, opacity, label, tipText, markers) {
  const exons = ivs(t.exons), cds = ivs(t.cds), mid = y + 9;
  let g = `<g opacity="${opacity}" data-tip="${esc(tipText)}"><line x1="${clip(t.start)}" y1="${mid}" x2="${clip(t.end)}" y2="${mid}" stroke="${color}" stroke-width="1"/>`;
  // strand chevrons
  const s0 = clip(t.start), s1 = clip(t.end);
  for (let px = s0 + 18; px < s1 - 10; px += 60) g += `<text x="${px}" y="${mid + 3.5}" font-size="9" fill="${color}" text-anchor="middle">${t.strand === "-" ? "‹" : "›"}</text>`;
  for (const [s, e] of exons) g += `<rect x="${clip(s)}" y="${mid - 4}" width="${Math.max(1, clip(e + 1) - clip(s))}" height="8" fill="var(--surface)" stroke="${color}" stroke-width="1.2"/>`;
  for (const [s, e] of cds) g += `<rect x="${clip(s)}" y="${mid - 7}" width="${Math.max(1, clip(e + 1) - clip(s))}" height="14" fill="${color}"/>`;
  g += `<text x="${4 + (markers.length ? 10 * markers.length + 2 : 0)}" y="${mid + 4}" class="mono" font-size="10">${esc(label.length > 38 ? label.slice(0, 37) + "…" : label)}</text>`;
  markers.forEach((c, i) => (g += `<rect x="${4 + i * 10}" y="${mid - 4}" width="8" height="8" fill="${c || "var(--na)"}" transform="rotate(45 ${8 + i * 10} ${mid})"/>`));
  return g + "</g>";
}
async function loadGeneDetail(geneId) {
  let d;
  try { d = await api("/api/gene", { reference: S.ref, policy: S.policy, gene_id: geneId }); } catch (e) { showError($("#geneDetail"), e); return; }
  const g = d.gene, runs = selectedRuns();
  const rows = runs.map((r) => {
    const o = d.outcomes[r.run_id]; const res = result(r.run_id, S.policy);
    if (!o) return `<tr><td>${runLabel(r)}</td><td colspan="7" class="na">${esc(res && res.status !== "complete" ? RESULT_LABEL[res.status] + (res.reason ? " — " + res.reason : "") : "no outcome")}</td></tr>`;
    const ex = o.explanation || {};
    const why = {
      exact: "Best CDS pair has identical CDS intervals.", terminal_diff: "Same CDS intron chain with ≥ 0.8 reciprocal overlap, but start and/or stop coordinates differ.",
      diff_chain: "≥ 0.8 reciprocal CDS overlap but the CDS intron chains differ (splice site or exon difference).",
      partial_lt08: "Shares coding sequence, but reciprocal CDS overlap is below 0.8.", span_only: "Gene spans overlap on the same strand without shared CDS.",
      strand_mismatch: `Only opposite-strand partners overlap (${o.strand_mismatch_basis}).`, missed: "No same-strand partner and no opposite-strand feature overlap.",
    }[o.cds_status] || "";
    return `<tr><td>${runLabel(r)}</td><td>${statusChip(o)}<div class="cell-sub">${esc(o.classification)} / CDS ${esc(o.classification_cds)}</div></td>
      <td class="mono small">${esc(o.best_cds_ref_tx || "—")}</td><td class="mono small">${esc(o.cds_matched_id || o.matched_id || "—")}<div class="cell-sub">${esc(o.best_cds_query_tx || "")}</div></td>
      <td class="small">${o.one_to_one_partner ? `<span class="status st-exact">${esc(o.one_to_one_partner)}</span>` : o.cds_exact ? '<span class="tag warn" data-tip="coordinate-exact, but the maximum one-to-one matching used this query gene (or this gene) in another pair">not in 1:1</span>' : '<span class="muted">—</span>'}${o.query_stop_codon ? `<div class="cell-sub" data-tip="${esc(`stop codon audit for ${o.query_stop_codon.transcript_id}: ${o.query_stop_codon.reason}\nlast codon ${o.query_stop_codon.last_codon || "?"}, next ${o.query_stop_codon.next_codon || "?"}`)}">stop ${esc(o.query_stop_codon.status)}: ${esc(o.query_stop_codon.terminal_cds_start)}-${esc(o.query_stop_codon.terminal_cds_end)}${o.query_stop_codon.status === "extended" ? " → " + esc(o.query_stop_codon.new_terminal_cds_start) + "-" + esc(o.query_stop_codon.new_terminal_cds_end) : ""}</div>` : ""}</td>
      <td class="num">${Number(o.cds_overlap).toFixed(3)}</td>
      <td class="small">best pair: ${esc(o.cds_intron_chain_match || "NA")}<br>any pair: ${o.any_pair_chain === null || o.any_pair_chain === undefined ? '<span class="na">n/a</span>' : o.any_pair_chain ? "yes" : "no"}${o.counterpart_count >= 2 ? `<br><span class="tag warn">split</span> ${esc(o.counterpart_ids)}` : ""}</td>
      <td class="small">${esc(why)}${ex.available ? `<div class="cell-sub">${esc(ex.summary)}${ex.reference_introns_missing_in_query?.length ? "<br>ref introns missing: " + esc(ex.reference_introns_missing_in_query.join(", ")) : ""}${ex.query_introns_not_in_reference?.length ? "<br>query-only introns: " + esc(ex.query_introns_not_in_reference.join(", ")) : ""}</div>` : ""}</td></tr>`;
  }).join("");
  $("#geneDetail").innerHTML = `<h3>${esc(g.gene_id)} ${g.name ? "· " + esc(g.name) : ""} <span class="small muted">${esc(g.biotype)} · ${esc(g.chrom)}:${fmtInt(g.start)}-${fmtInt(g.end)} ${esc(g.strand)}${g.original_gene_id !== g.gene_id ? " · file ID " + esc(g.original_gene_id) : ""}</span></h3>
    <p class="small muted">Policy ${esc(S.policy)}: ${g.policy ? `${g.policy.evaluated_transcripts.split(",").length} evaluated isoform(s): ${esc(g.policy.evaluated_transcripts.split(",").slice(0, 8).join(", "))}${g.policy.evaluated_transcripts.split(",").length > 8 ? ", … (" + (g.policy.evaluated_transcripts.split(",").length - 8) + " more)" : ""}${g.policy.canonical_fallback ? " (no canonical tag: longest CDS)" : ""}; ${g.policy.multi_segment ? "multi" : "single"}-segment CDS` : "not evaluated under this policy"}. Differences are recomputed from the stored models (positive bp = query CDS extends further).</p>
    <table class="grid"><tr><th>Run</th><th>Outcome</th><th>Compared reference isoform</th><th>Query partner</th><th>One-to-one exact pair</th><th>CDS reciprocal overlap</th><th>Intron chain</th><th>Why</th></tr>${rows}</table>`;
}

// ---------------------------------------------------------------- provenance
function kv(obj) { return `<dl class="kv">${Object.entries(obj).map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v === null || v === undefined || v === "" ? '<span class="na">unknown</span>' : esc(typeof v === "object" ? JSON.stringify(v) : v)}</dd>`).join("")}</dl>`; }
function renderProvenance() {
  const m = S.meta, el = $("#provenance"), v = m.validation;
  const ref = m.references.find((r) => r.id === S.ref), prep = ref.preparation;
  const issues = v ? v.issues.filter((i) => i.where.includes(S.ref) || refRuns().some((r) => i.where.includes(r.run_id))) : [];
  let html = `<h2>Experiment</h2>${kv({ id: m.experiment.id, title: m.experiment.title, synthetic: m.experiment.synthetic ? "YES — synthetic demonstration data" : "no",
    configuration: m.experiment.config_path, "dataset built": m.built, "validated": v?.validated || "not validated" })}`;
  html += `<h2>Validation for this reference ${v ? (v.ok ? '<span class="tag good">no errors</span>' : '<span class="tag bad">errors</span>') : '<span class="tag warn">not run</span>'}</h2>`;
  html += issues.length ? issues.map((i) => `<div class="issue ${esc(i.level)}"><b>${esc(i.level)}</b> · ${esc(i.where)}: ${esc(i.message)}</div>`).join("") : '<p class="muted">No issues recorded for this reference and its runs.</p>';
  if (v) {
    html += `<h3>Known comparator parser limitations (probed against the checkout at validation time)</h3><table class="grid"><tr><th>Limitation</th><th>Status</th><th>Handling</th></tr>${v.parser_limitations.map((l) => `<tr><td>${esc(l.title)}</td><td>${l.status.startsWith("present") ? '<span class="tag warn">' : '<span class="tag good">'}${esc(l.status)}</span></td><td class="small">${esc(l.handling)}</td></tr>`).join("")}</table>`;
  }
  html += `<h2>Reference ${esc(ref.id)}</h2>${kv({ species: ref.species, assembly: ref.assembly, accession: ref.assembly_accession, source: ref.annotation_source, release: ref.annotation_release, notes: ref.notes })}`;
  if (prep) {
    for (const k of ["annotation", "genome"]) {
      const p = prep[k]; if (!p) continue;
      html += `<h3>${k}: ${esc(p.method)}</h3>${kv({ source: p.source.path, "source sha256": p.source.sha256, recipe: p.recipe.map((s) => s.name + (s.types ? "(" + s.types.join(",") + ")" : "")).join(" → ") || "none", output: p.output.path, "output sha256": p.output.sha256, counts: p.counts && Object.keys(p.counts).length ? p.counts : "none" })}`;
    }
    const sv = v?.references?.[S.ref];
    html += `<h3>Scope</h3>${kv({ type: prep.scope.type, sequences: prep.scope.sequences ? prep.scope.sequences.length : "all", "regions file": prep.scope.regions_file, "sha256": prep.scope.sha256,
      "canonical tags present": sv ? String(sv.canonical_available) : "", "genome sequences": sv?.genome_sequences })}`;
  } else html += '<p class="na">Reference not prepared in this workspace: locus browsing and segment-level metrics unavailable.</p>';
  html += "<h2>Prediction runs</h2>";
  for (const r of refRuns()) {
    const rv = v?.runs?.[r.run_id] || {};
    const ms = r.models_source || {};
    html += `<h3>${runLabel(r)}</h3>${kv({ label: runName(r), "run id": r.run_id, tool: r.tool, role: r.role, "listed after": r.companion_of ? runDisplay(runById(r.companion_of)) : "", "query coordinates": harmonised(r) ? "stop-harmonised (--add-stop-codon)" : "original",
      "comparator options": Object.keys(r.comparator_options || {}).length ? r.comparator_options : "none (original coordinates)",
      "locus models from": ms.kind ? `${ms.kind.replaceAll("_", " ")}${ms.path ? ": " + ms.path : ""}${ms.policy_id ? " (policy " + ms.policy_id + " output)" : ""}${ms.reason ? " — " + ms.reason : ""}` : "",
      version: r.tool_version, model: r.model, annotation: r.annotation, "input sha256": rv.input?.sha256, format: rv.annotation_scan?.format,
      provenance: r.provenance, "predictions outside scope": rv.outside_scope ? `${fmtInt(rv.outside_scope.genes)} genes on ${rv.outside_scope.sequences} sequences (excluded)` : "",
      "scope sequences without predictions": rv.scope_sequences_without_predictions, "locus models": m.model_errors?.[r.run_id] ? "unavailable: " + m.model_errors[r.run_id] : "available" })}`;
    if ((r.caveats || []).length) html += `<div class="issue warning"><b>Caveats</b><ul>${r.caveats.map((c) => `<li>${esc(c)}</li>`).join("")}</ul></div>`;
    const stop = m.policies.map((p) => result(r.run_id, p.id)?.stop_codon_harmonisation).find((b) => b && b.applied);
    if (stop) html += `<div class="issue warning"><b>Stop-codon harmonisation (query)</b> · ${fmtInt(stop.transcripts_with_cds)} transcripts with CDS: ${fmtInt(stop.eligible)} frame established, ${fmtInt(stop.extended)} extended, ${fmtInt(stop.unchanged)} unchanged, ${fmtInt(stop.skipped)} skipped
      <ul>${Object.entries(stop.reasons || {}).map(([k, v]) => `<li class="small">${esc(k)}: ${fmtInt(v)}</li>`).join("")}</ul><span class="small">Genetic code: default table ${esc(stop.policy?.genetic_code?.default_table)}; overrides ${esc(JSON.stringify(stop.policy?.genetic_code?.sequence_overrides || {}))}; organelle names need an explicit override. Phase not used.</span></div>`;
    html += `<table class="grid"><tr><th>Policy</th><th>Status</th><th>Source</th><th>Independent</th><th>1:1 pairs export</th><th>Finished</th><th>Wall s</th><th>Peak RSS</th><th>Output</th></tr>${m.policies.map((p) => {
      const res = result(r.run_id, p.id) || {};
      return `<tr><td>${esc(p.id)}</td><td>${esc(RESULT_LABEL[res.status] || res.status || "no result")}${res.reason && res.status !== "complete" ? `<div class="cell-sub">${esc(res.reason)}</div>` : ""}${res.difference ? `<div class="cell-sub">${esc(res.difference)}</div>` : ""}</td>
        <td>${esc(res.source || "")}</td><td class="small">${esc(res.independent?.status || "")}${res.independent?.reason ? " — " + esc(res.independent.reason) : ""}</td>
        <td class="small">${res.pairs_check ? `${esc(res.pairs_check.status)}${res.pairs_check.rows !== undefined ? ` (${fmtInt(res.pairs_check.rows)} rows, TP ${fmtInt(res.pairs_check.true_positives)})` : ""}` : ""}</td><td class="small">${esc(res.finished || "")}</td>
        <td class="num">${res.wall_seconds ?? ""}</td><td class="num">${res.peak_rss_bytes ? (res.peak_rss_bytes / 1e9).toFixed(2) + " GB" : '<span class="na">not recorded</span>'}</td><td class="mono small">${esc(res.output_dir || "")}</td></tr>`;
    }).join("")}</table>`;
  }
  const imp = m.import_report;
  if (imp) {
    html += `<h2>Import</h2>${kv({ "benchmark directory": imp.benchmark_dir, imported: imp.imported, "whole-genome results matched": imp.matched.length, "pilots (kept separate)": imp.pilots.length, unmatched: imp.unmatched.length, conflicts: imp.conflicts.length, "code snapshot": imp.code_snapshot.join(", ") })}`;
    if (imp.pilots.length) html += `<details><summary>Region pilots — not whole-genome results, not shown in the overview</summary><table class="grid"><tr><th>Directory</th><th>Run</th><th>Region</th><th>Reason</th></tr>${imp.pilots.map((p) => `<tr><td class="mono small">${esc(p.dir)}</td><td>${esc(p.run_id)}</td><td>${esc((p.options_region || []).join(","))}</td><td class="small">${esc(p.reason)}</td></tr>`).join("")}</table></details>`;
    if (imp.unmatched.length) html += `<details><summary>Unmatched directories</summary><pre class="mono">${esc(JSON.stringify(imp.unmatched, null, 1))}</pre></details>`;
  }
  el.innerHTML = html;
}

// ---------------------------------------------------------------- dictionary & context
function renderDictionary() {
  const d = S.meta.metric_dictionary;
  $("#dictionary").innerHTML = `<h2>Metric dictionary</h2><ul>${d.notes.map((n) => `<li class="small">${esc(n)}</li>`).join("")}</ul>
    <table class="grid"><tr><th>Metric</th><th>Group</th><th>Unit</th><th>Numerator</th><th>Denominator</th><th>Source</th><th>Source fields</th><th>Limitations</th></tr>
    ${d.metrics.map((m) => `<tr><td><b>${esc(m.label)}</b><div class="cell-sub mono">${esc(m.id)}</div>${m.definition ? `<div class="small">${esc(m.definition)}</div>` : ""}</td><td>${esc(m.group)}</td><td>${esc(m.unit)}</td>
      <td class="small">${esc(m.numerator)}</td><td class="small">${esc(m.denominator)}</td><td><span class="tag">${esc(m.source)}</span></td>
      <td class="small mono">${esc((m.source_fields || []).join("; "))}</td><td class="small">${(m.limitations || []).map((l) => "• " + esc(l)).join("<br>")}</td></tr>`).join("")}</table>
    <h3>${esc(d.paper_analogous.label)}</h3><p class="small">${esc(d.paper_analogous.definition)}</p><ul>${d.paper_analogous.limitations.map((l) => `<li class="small">${esc(l)}</li>`).join("")}</ul>`;
}
function renderContext() {
  const m = S.meta;
  let html = `<div class="callout"><b>Paper-analogous metrics are supporting context, not a reproduction</b> of the Tiberius benchmark (different reference, tool, matching rules and model versions; mouse was a Tiberius training genome). Verified methods and stated assumptions are distinguished in the method-comparison document below.</div>`;
  for (const p of m.paper_analogous || []) {
    if (!p.data) { html += `<p class="na">${esc(p.name)} not available.</p>`; continue; }
    const policyLabel = /refB/.test(p.name) ? "reference = longest CDS per gene (the paper's policy)" : /refC/.test(p.name) ? "reference = Ensembl_canonical (supplementary)" : "";
    html += `<h3>Paper-analogous · ${esc(policyLabel)} <span class="small muted">${esc(p.name)} · query = longest CDS · one transcript per gene · one-to-one matching</span></h3><table class="grid"><tr><th>Run (file label)</th><th>Gene F1 coordinate-exact</th><th>Gene F1 intron-chain (gffcompare-like)</th><th>CDS-exon F1</th><th>R / Q</th></tr>${Object.entries(p.data).map(([k, v]) =>
      `<tr><td>${esc(k)}</td><td class="num">${v.gene_exact.F1.toFixed(3)}</td><td class="num">${v.gene_chain.F1.toFixed(3)}</td><td class="num">${v.exon_cds_segments.F1.toFixed(3)}</td><td class="num">${fmtInt(v.reference_genes)} / ${fmtInt(v.query_genes)}</td></tr>`).join("")}</table>`;
  }
  html += `<h2>Documents</h2>${(m.context_documents || []).map((d) => `<button class="btn" data-doc="${d.index}" ${d.available ? "" : "disabled"}>${esc(d.name)}</button>`).join(" ") || '<p class="muted">No context documents configured.</p>'}<div id="docView" class="md"></div>`;
  $("#context").innerHTML = html;
  $$("[data-doc]").forEach((b) => (b.onclick = async () => { const d = await api("/api/document", { index: b.dataset.doc }); $("#docView").innerHTML = `<p class="small muted">${esc(d.path)}</p>` + markdown(d.text || ""); }));
}
function markdown(text) {
  const lines = text.split("\n"); let out = "", i = 0, list = null;
  const inline = (s) => esc(s).replace(/`([^`]+)`/g, "<code>$1</code>").replace(/\*\*([^*]+)\*\*/g, "<b>$1</b>").replace(/\*([^*]+)\*/g, "<i>$1</i>")
    .replace(/!\[([^\]]*)\]\([^)]*\)/g, "[image: $1]").replace(/\[([^\]]+)\]\((https?:[^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>').replace(/\[([^\]]+)\]\(([^)]+)\)/g, "$1");
  const closeList = () => { if (list) { out += `</${list}>`; list = null; } };
  while (i < lines.length) {
    const l = lines[i];
    if (l.startsWith("```")) { closeList(); let code = ""; i++; while (i < lines.length && !lines[i].startsWith("```")) code += lines[i++] + "\n"; out += `<pre>${esc(code)}</pre>`; i++; continue; }
    if (/^\|.*\|\s*$/.test(l)) { closeList(); const rows = []; while (i < lines.length && /^\|.*\|\s*$/.test(lines[i])) rows.push(lines[i++]);
      out += "<table>" + rows.filter((r) => !/^\|\s*:?-+/.test(r)).map((r, k) => "<tr>" + r.trim().slice(1, -1).split("|").map((c) => (k === 0 ? `<th>${inline(c.trim())}</th>` : `<td>${inline(c.trim())}</td>`)).join("") + "</tr>").join("") + "</table>"; continue; }
    const h = /^(#{1,4})\s+(.*)/.exec(l);
    if (h) { closeList(); out += `<h${h[1].length + 1}>${inline(h[2])}</h${h[1].length + 1}>`; i++; continue; }
    const li = /^\s*(?:[-*]|\d+\.)\s+(.*)/.exec(l);
    if (li) { const type = /^\s*\d+\./.test(l) ? "ol" : "ul"; if (list !== type) { closeList(); out += `<${type}>`; list = type; } out += `<li>${inline(li[1])}</li>`; i++; continue; }
    closeList();
    if (l.trim()) out += `<p>${inline(l)}</p>`;
    i++;
  }
  closeList();
  return out;
}

// Deep links: #exp=…&ref=…&policy=… (applied on load and when the hash changes)
window.addEventListener("hashchange", () => {
  const p = new URLSearchParams(location.hash.slice(1));
  const exp = p.get("exp"), ref = p.get("ref"), policy = p.get("policy"), scope = p.get("scope");
  if (exp && exp !== S.exp) { S.exp = exp; $("#exp").value = exp; loadExperiment(p); return; }
  let changed = false;
  if (policy && policy !== S.policy && S.meta.policies.some((x) => x.id === policy)) { S.policy = policy; changed = true; }
  if (ref && ref !== S.ref && S.meta.references.some((r) => r.id === ref)) { S.ref = ref; $("#ref").value = ref; S.runs.clear(); loadReference(); return; }
  if (scope && scope !== S.scope && scopes().some((x) => x.id === scope && x.available)) { S.scope = scope; changed = true; }
  if (changed) { renderPolicySelector(); renderScopeSelector(); renderRuns(); renderOverview(); S.genesLoaded = S.predsLoaded = false; }
});

init().catch((e) => { document.body.insertAdjacentHTML("afterbegin", `<div class="issue error">${esc(e.message)}</div>`); });
