// Stage colours of the reaction workflow (stages and pages: p2r_pipeline in conf.py).
// With html.p2r-deco (conf.py p2r_deco), adds a stage dot to command links in the
// sidebar and in command tables, and colours the numbered "How it works" steps on
// the `all` page.
(function () {
  "use strict";

  if (!document.documentElement.classList.contains("p2r-deco")) return;

  var stageOf = window.p2rStageOf || {};
  var order = window.p2rStageOrder || [];

  function rootUrl() {
    var r = document.documentElement.getAttribute("data-content_root") || "./";
    return new URL(r, location.href).href;
  }

  // Page name of a link relative to the docs root, without "ja/" and ".html".
  function pageOf(a, root) {
    var href = a.getAttribute("href");
    if (href === null) return null;
    var url = new URL(href, location.href).href.split("#")[0];
    if (url.indexOf(root) !== 0) return null;
    var rel = url.slice(root.length).replace(/\.html$/, "");
    return rel.indexOf("ja/") === 0 ? rel.slice(3) : rel;
  }

  function mark(el, stage) {
    el.setAttribute("data-p2r-stage", stage);
    if (order.indexOf(stage) >= 0) el.style.setProperty("--p2r-c", "var(--p2r-" + stage + ")");
  }

  // Mark every link of a list once any of them is a command with a stage, so the
  // neutral commands of the same list get an empty dot and the column stays aligned.
  function markGroup(links, root) {
    var stages = links.map(function (a) { return stageOf[pageOf(a, root)] || null; });
    if (!stages.some(Boolean)) return;
    links.forEach(function (a, i) { mark(a, stages[i] || "none"); });
  }

  function run() {
    var root = rootUrl();

    document.querySelectorAll(".sidebar-tree ul").forEach(function (ul) {
      var links = [];
      ul.querySelectorAll(":scope > li.toctree-l1 > a.reference").forEach(function (a) { links.push(a); });
      markGroup(links, root);
    });

    document.querySelectorAll("article table").forEach(function (table) {
      var links = [];
      table.querySelectorAll("tbody > tr > td:first-child").forEach(function (td) {
        var a = td.querySelector("a.reference.internal");
        if (a) links.push(a);
      });
      markGroup(links, root);
    });

    var all = document.querySelector('article[data-p2r-stage="all"]');
    if (!all) return;
    var steps = window.p2rAllSteps || [];
    var lists = all.querySelectorAll("section > ol");
    for (var k = 0; k < lists.length; k++) {
      var items = lists[k].querySelectorAll(":scope > li");
      if (!steps.length || items.length !== steps.length) continue;
      var ok = true;
      items.forEach(function (li) {
        if (!li.querySelector(":scope > p:first-child > strong:first-child")) ok = false;
      });
      if (!ok) continue;
      lists[k].classList.add("p2r-steps");
      // A step that covers two stages (thermochemistry and DFT) gets a two-colour disc.
      items.forEach(function (li, i) {
        var ids = steps[i];
        li.style.setProperty("--p2r-c", "var(--p2r-" + ids[0] + ")");
        if (ids.length > 1) {
          li.style.setProperty("--p2r-disc", "linear-gradient(135deg, " + ids.map(function (id, j) {
            return "var(--p2r-" + id + ") " + (100 * j / ids.length) + "% " + (100 * (j + 1) / ids.length) + "%";
          }).join(", ") + ")");
        }
      });
      break;
    }
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", run);
  } else {
    run();
  }
})();
