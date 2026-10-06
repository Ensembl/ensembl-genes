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
const STRUCTURE = ["cds_exact_ref_multi", "cds_exact_ref_single", "partial_cds", "partial_same_chain_terminal_diff", "partial_ge08_diff_chain",
  "partial_lt08", "cds_structural_exact_match", "gene_span_locus_detection", "gene_span_only_no_exon_overlap", "detected_without_cds_overlap",
  "query_matched_without_cds_overlap", "exact_tp_one_to_one", "exact_pairs_query_gene_reused"];
const LOCI = ["missed_reference_genes", "ref_genes_without_cds_overlap", "novel_predictions", "splits", "merges", "one_to_one_ref_genes",
  "strand_mismatch_ref", "strand_mismatch_query"];
const NOVEL = ["novel_Overlaps_reference_pseudogene", "novel_Overlaps_other_non_protein_coding_reference",
  "novel_Overlaps_unevaluated_protein_coding_reference", "novel_Overlaps_evaluated_reference_opposite_strand", "novel_Novel_no_reference_overlap"];
const NOVEL_COLORS = ["#184f95", "#2a78d6", "#6da7ec", "#9ec5f4", "#cde2fb"];
const UTR = ["exon_coordinate_exact_ref", "structural_exact_match_with_utr"];
const SCALE = ["reference_coding_genes", "reference_eligible_transcripts", "predicted_coding_genes", "predicted_transcripts"];

const S = { exp: null, meta: null, ref: null, policy: null, runs: new Set(), ov: null, dict: {}, genePage: 0, predPage: 0,
  locus: null, selectedGene: null, pageSize: 100 };

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
  await loadExperiment(params);
}
function selectTab(name) {
  $$(".tabs button").forEach((b) => b.classList.toggle("active", b.dataset.tab === name));
  $$(".tab").forEach((t) => t.classList.toggle("active", t.id === "tab-" + name));
  if (name === "genes" && !S.genesLoaded) loadGenes();
  if (name === "predictions" && !S.predsLoaded) loadPredictions();
}
async function loadExperiment(params = new URLSearchParams()) {
  S.meta = await api("/api/meta");
  S.dict = Object.fromEntries((S.meta.metric_dictionary?.metrics || []).map((m) => [m.id, m]));
  $("#synthetic").hidden = !S.meta.experiment.synthetic;
  const refs = S.meta.references;
  $("#ref").innerHTML = refs.map((r) => `<option value="${esc(r.id)}">${esc(r.species)} · ${esc(r.assembly)} · ${esc(r.id)}</option>`).join("");
  S.ref = params.get("ref") && refs.some((r) => r.id === params.get("ref")) ? params.get("ref") : refs[0].id;
  $("#ref").value = S.ref;
  $("#ref").onchange = () => { S.ref = $("#ref").value; S.runs.clear(); loadReference(); };
  S.policy = params.get("policy") || S.meta.policies[0].id;
  renderDictionary(); renderContext();
  await loadReference();
}
function refRuns() { return S.meta.runs.filter((r) => r.reference === S.ref); }
function runById(id) { return S.meta.runs.find((r) => r.run_id === id); }
function selectedRuns() { return refRuns().filter((r) => S.runs.has(r.run_id)); }
function result(run, policy) { return S.ov.results.find((r) => r.run_id === run && r.policy_id === policy); }
function metric(run, policy, id) { return S.ov.metrics.find((m) => m.run_id === run && m.policy_id === policy && m.metric_id === id); }
function runLabel(r, html = true) {
  const detail = [r.model ? `model ${r.model}` : "model unknown", r.tool_version ? `v${r.tool_version}` : "version unknown"].join(", ");
  if (!html) return `${r.tool} (${r.run_id})`;
  return `<span class="swatch" style="background:${r.color}"></span><span><b>${esc(r.tool)}</b> <span class="small muted" data-tip="${esc(r.run_id + "\n" + detail)}">${esc(r.run_id)}</span></span>`;
}

async function loadReference() {
  S.ov = await api("/api/overview", { reference: S.ref });
  if (!S.runs.size) refRuns().forEach((r) => S.runs.add(r.run_id));
  renderPolicySelector(); renderRuns(); renderOverview(); renderProvenance();
  S.genesLoaded = S.predsLoaded = false;
  $("#locusPanel").hidden = true;
  setupGeneFilters(); setupPredFilters();
  const active = $(".tabs button.active").dataset.tab;
  if (active === "genes") loadGenes();
  if (active === "predictions") loadPredictions();
  writeHash();
}
function writeHash() {
  history.replaceState(null, "", `#exp=${encodeURIComponent(S.exp)}&ref=${encodeURIComponent(S.ref)}&policy=${encodeURIComponent(S.policy)}`);
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
      data-tip="${esc(p.label + "\n" + (p.description || "") + (av.ok ? "" : "\nUnavailable: " + av.why))}">${esc(p.id)} · ${esc(p.label)}</button>`;
  }).join("");
  $$("button", box).forEach((b) => (b.onclick = () => { S.policy = b.dataset.p; renderPolicySelector(); renderOverview(); S.genesLoaded = S.predsLoaded = false;
    const t = $(".tabs button.active").dataset.tab; if (t === "genes") loadGenes(); if (t === "predictions") loadPredictions(); if (S.locus) drawLocus(); writeHash(); }));
  const p = S.meta.policies.find((x) => x.id === S.policy);
  const av = policyAvailability(S.policy);
  $("#policyNote").innerHTML = `${esc(p ? p.description : "")}${av.ok ? "" : ` <span class="tag warn">unavailable for this reference: ${esc(av.why)}</span>`}`;
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
  if (!m) return `<td class="num"><span class="na">—</span></td>`;
  if (m.status !== "ok") return `<td class="num"><span class="na" data-tip="${esc(m.note || "")}">n/a</span></td>`;
  const def = S.dict[id] || {};
  if (def.unit === "rate") return `<td class="num"><div class="cell-main">${(m.value ?? 0).toFixed(3)}</div><div class="cell-sub">TP ${fmtInt(m.numerator)}${m.denominator ? " / " + fmtInt(m.denominator) : ""}</div></td>`;
  if (m.denominator === null || m.denominator === undefined) return `<td class="num"><div class="cell-main">${fmtInt(m.numerator)}</div></td>`;
  return `<td class="num"><div class="cell-main">${fmtPct(m.numerator / m.denominator)}</div><div class="cell-sub">${fmtInt(m.numerator)} / ${fmtInt(m.denominator)}</div></td>`;
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
    rows.push([id, S.dict[id]?.label || id, S.dict[id]?.source || "", r.run_id, r.tool, S.policy, res?.status || "no_result", m?.status || "",
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
  const ref = S.meta.references.find((r) => r.id === S.ref);
  const prep = ref.preparation || {};
  const R = refRuns().map((r) => metric(r.run_id, S.policy, "reference_coding_genes")).find((m) => m && m.status === "ok");
  const scope = prep.scope ? (prep.scope.type === "all" ? "all sequences" : `${(prep.scope.sequences || []).length} sequences (${prep.scope.type})`) : "unknown (not prepared)";
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
  const sel = $("#policyMetric");
  const rateIds = HEADLINE.concat(STRUCTURE).filter((id) => S.dict[id] && S.dict[id].denominator && S.dict[id].denominator !== "—");
  if (!sel.options.length) { sel.innerHTML = rateIds.map((id) => `<option value="${id}">${esc(S.dict[id].label)}</option>`).join(""); sel.value = "cds_exact_ref"; sel.onchange = renderPolicyChart; }
  renderPolicyChart(); renderPartialChart(); renderNovelChart();
  $("#dlMetrics").onclick = () => (location.href = apiUrl("/api/overview", { reference: S.ref, format: "csv" }));
}

// ---------------------------------------------------------------- charts (SVG)
function svg(w, h, inner) { return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img">${inner}</svg>`; }
function renderPolicyChart() {
  const id = $("#policyMetric").value; const runs = selectedRuns(); const pols = S.meta.policies;
  const el = $("#policyChart");
  if (!runs.length) { el.innerHTML = ""; return; }
  const W = Math.min(980, el.clientWidth || 980), L = 210, bh = 16, gap = 4, groupGap = 14;
  let y = 8, out = "";
  const x = (v) => L + v * (W - L - 130);
  for (const p of pols) {
    out += `<text x="4" y="${y + 12}" font-weight="600">${esc(p.id)} · ${esc(p.label)}</text>`;
    y += 18;
    for (const r of runs) {
      const res = result(r.run_id, p.id); const m = metric(r.run_id, p.id, id);
      out += `<text x="14" y="${y + 12}" class="lbl-muted">${esc(r.tool)} <tspan fill="var(--muted)">${esc(r.run_id)}</tspan></text>`;
      if (!res || res.status !== "complete") {
        out += `<text x="${L}" y="${y + 12}" fill="var(--na)" font-style="italic">${esc(res ? RESULT_LABEL[res.status] : "no result")}${res?.reason ? " — " + esc(res.reason.slice(0, 80)) : ""}</text>`;
      } else if (!m || m.status !== "ok") {
        out += `<text x="${L}" y="${y + 12}" fill="var(--na)" font-style="italic">n/a — ${esc(m?.note || "")}</text>`;
      } else {
        const v = m.numerator / m.denominator;
        out += `<rect x="${L}" y="${y + 2}" width="${Math.max(1, x(v) - L)}" height="${bh - 4}" rx="3" fill="${r.color}" data-tip="${esc(`${runLabel(r, false)}, policy ${p.id}\n${fmtInt(m.numerator)} / ${fmtInt(m.denominator)} = ${fmtPct(v, 2)}`)}"></rect>`;
        out += `<text x="${x(v) + 6}" y="${y + 12}">${fmtPct(v)} <tspan fill="var(--muted)">(${fmtInt(m.numerator)} / ${fmtInt(m.denominator)})</tspan></text>`;
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
  const W = Math.min(980, el.clientWidth || 980), L = 210, bh = 18;
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
    out += `<text x="4" y="${y + 13}">${esc(r.tool)} <tspan fill="var(--muted)">${esc(r.run_id)}</tspan></text>`;
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
  $("#focusRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(r.tool)} · ${esc(r.run_id)}</option>`).join("");
  $("#compareRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(r.tool)} · ${esc(r.run_id)}</option>`).join("");
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
    o.cds_matched_id ? `partner ${o.cds_matched_id}` : "", o.counterpart_count >= 2 ? `split: ${o.counterpart_ids}` : "",
    o.any_pair_chain === null || o.any_pair_chain === undefined ? "" : `any-pair intron chain (independent): ${o.any_pair_chain ? "yes" : "no"}`].filter(Boolean).join("\n");
  return `<span class="status st-${esc(o.cds_status)}" data-tip="${esc(t + "\n" + extra)}">${esc(lab)}${o.counterpart_count >= 2 ? " ⇉" : ""}</span>`;
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
  $("#predRun").innerHTML = runs.map((r) => `<option value="${esc(r.run_id)}">${esc(r.tool)} · ${esc(r.run_id)}</option>`).join("");
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
    out += `<rect x="4" y="${y + 3}" width="10" height="10" fill="${r.color}" rx="2"/><text x="18" y="${y + 12}" font-weight="600">${esc(r.tool)} <tspan class="lbl-muted" font-weight="400">${esc(r.run_id)}</tspan></text>`;
    y += 18;
    const modelError = S.meta.model_errors?.[r.run_id];
    if (modelError || !tr || !tr.transcripts.length) {
      const msg = modelError ? `prediction models unavailable — ${modelError.slice(0, 110)}` : tr ? "no predictions in this window" : "prediction models unavailable";
      out += `<text x="${LBL}" y="${y + 12}" class="lbl-muted">${esc(msg)}</text>`; y += rowH; continue;
    }
    for (const t of tr.transcripts) {
      const o = t.outcome;
      const lab = `${t.gene_id}${o ? " · " + o.classification + (o.cds_exact ? " · CDS exact" : "") : ""}`;
      out += transcriptGlyph(t, y, x, clip, r.color, 1, lab, `${t.transcript_id}\ngene ${t.gene_id}\n${o ? `comparator: ${o.classification} / CDS ${o.classification_cds}${o.novel_category ? "\nNovel: " + o.novel_category : ""}${o.counterpart_count >= 2 ? "\nmerge (" + o.counterpart_count + " reference counterparts)" : ""}` : "no outcome for this policy"}`, []);
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
    if (!o) return `<tr><td>${runLabel(r)}</td><td colspan="6" class="na">${esc(res && res.status !== "complete" ? RESULT_LABEL[res.status] + (res.reason ? " — " + res.reason : "") : "no outcome")}</td></tr>`;
    const ex = o.explanation || {};
    const why = {
      exact: "Best CDS pair has identical CDS intervals.", terminal_diff: "Same CDS intron chain with ≥ 0.8 reciprocal overlap, but start and/or stop coordinates differ.",
      diff_chain: "≥ 0.8 reciprocal CDS overlap but the CDS intron chains differ (splice site or exon difference).",
      partial_lt08: "Shares coding sequence, but reciprocal CDS overlap is below 0.8.", span_only: "Gene spans overlap on the same strand without shared CDS.",
      strand_mismatch: `Only opposite-strand partners overlap (${o.strand_mismatch_basis}).`, missed: "No same-strand partner and no opposite-strand feature overlap.",
    }[o.cds_status] || "";
    return `<tr><td>${runLabel(r)}</td><td>${statusChip(o)}<div class="cell-sub">${esc(o.classification)} / CDS ${esc(o.classification_cds)}</div></td>
      <td class="mono small">${esc(o.best_cds_ref_tx || "—")}</td><td class="mono small">${esc(o.cds_matched_id || o.matched_id || "—")}<div class="cell-sub">${esc(o.best_cds_query_tx || "")}</div></td>
      <td class="num">${Number(o.cds_overlap).toFixed(3)}</td>
      <td class="small">best pair: ${esc(o.cds_intron_chain_match || "NA")}<br>any pair: ${o.any_pair_chain === null || o.any_pair_chain === undefined ? '<span class="na">n/a</span>' : o.any_pair_chain ? "yes" : "no"}${o.counterpart_count >= 2 ? `<br><span class="tag warn">split</span> ${esc(o.counterpart_ids)}` : ""}</td>
      <td class="small">${esc(why)}${ex.available ? `<div class="cell-sub">${esc(ex.summary)}${ex.reference_introns_missing_in_query?.length ? "<br>ref introns missing: " + esc(ex.reference_introns_missing_in_query.join(", ")) : ""}${ex.query_introns_not_in_reference?.length ? "<br>query-only introns: " + esc(ex.query_introns_not_in_reference.join(", ")) : ""}</div>` : ""}</td></tr>`;
  }).join("");
  $("#geneDetail").innerHTML = `<h3>${esc(g.gene_id)} ${g.name ? "· " + esc(g.name) : ""} <span class="small muted">${esc(g.biotype)} · ${esc(g.chrom)}:${fmtInt(g.start)}-${fmtInt(g.end)} ${esc(g.strand)}${g.original_gene_id !== g.gene_id ? " · file ID " + esc(g.original_gene_id) : ""}</span></h3>
    <p class="small muted">Policy ${esc(S.policy)}: ${g.policy ? `${g.policy.evaluated_transcripts.split(",").length} evaluated isoform(s): ${esc(g.policy.evaluated_transcripts.split(",").slice(0, 8).join(", "))}${g.policy.evaluated_transcripts.split(",").length > 8 ? ", … (" + (g.policy.evaluated_transcripts.split(",").length - 8) + " more)" : ""}${g.policy.canonical_fallback ? " (no canonical tag: longest CDS)" : ""}; ${g.policy.multi_segment ? "multi" : "single"}-segment CDS` : "not evaluated under this policy"}. Differences are recomputed from the stored models (positive bp = query CDS extends further).</p>
    <table class="grid"><tr><th>Run</th><th>Outcome</th><th>Compared reference isoform</th><th>Query partner</th><th>CDS reciprocal overlap</th><th>Intron chain</th><th>Why</th></tr>${rows}</table>`;
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
    html += `<h3>${runLabel(r)}</h3>${kv({ tool: r.tool, version: r.tool_version, model: r.model, annotation: r.annotation, "input sha256": rv.input?.sha256, format: rv.annotation_scan?.format,
      provenance: r.provenance, "predictions outside scope": rv.outside_scope ? `${fmtInt(rv.outside_scope.genes)} genes on ${rv.outside_scope.sequences} sequences (excluded)` : "",
      "scope sequences without predictions": rv.scope_sequences_without_predictions, "locus models": m.model_errors?.[r.run_id] ? "unavailable: " + m.model_errors[r.run_id] : "available" })}`;
    if ((r.caveats || []).length) html += `<div class="issue warning"><b>Caveats</b><ul>${r.caveats.map((c) => `<li>${esc(c)}</li>`).join("")}</ul></div>`;
    html += `<table class="grid"><tr><th>Policy</th><th>Status</th><th>Source</th><th>Independent</th><th>Finished</th><th>Wall s</th><th>Peak RSS</th><th>Output</th></tr>${m.policies.map((p) => {
      const res = result(r.run_id, p.id) || {};
      return `<tr><td>${esc(p.id)}</td><td>${esc(RESULT_LABEL[res.status] || res.status || "no result")}${res.reason && res.status !== "complete" ? `<div class="cell-sub">${esc(res.reason)}</div>` : ""}${res.difference ? `<div class="cell-sub">${esc(res.difference)}</div>` : ""}</td>
        <td>${esc(res.source || "")}</td><td class="small">${esc(res.independent?.status || "")}${res.independent?.reason ? " — " + esc(res.independent.reason) : ""}</td><td class="small">${esc(res.finished || "")}</td>
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
  const exp = p.get("exp"), ref = p.get("ref"), policy = p.get("policy");
  if (exp && exp !== S.exp) { S.exp = exp; $("#exp").value = exp; loadExperiment(p); return; }
  let changed = false;
  if (policy && policy !== S.policy && S.meta.policies.some((x) => x.id === policy)) { S.policy = policy; changed = true; }
  if (ref && ref !== S.ref && S.meta.references.some((r) => r.id === ref)) { S.ref = ref; $("#ref").value = ref; S.runs.clear(); loadReference(); return; }
  if (changed) { renderPolicySelector(); renderOverview(); S.genesLoaded = S.predsLoaded = false; }
});

init().catch((e) => { document.body.insertAdjacentHTML("afterbegin", `<div class="issue error">${esc(e.message)}</div>`); });
