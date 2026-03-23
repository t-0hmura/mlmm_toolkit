// Override Furo's 3-state theme toggle (light/dark/auto) with 2-state (light ↔ dark).
(function () {
  "use strict";

  function resolveAuto() {
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }

  function applyTheme(theme) {
    document.body.dataset.theme = theme;
    localStorage.setItem("theme", theme);
  }

  function toggle() {
    var current = document.body.dataset.theme;
    if (current !== "light" && current !== "dark") current = resolveAuto();
    applyTheme(current === "dark" ? "light" : "dark");
  }

  // Intercept all theme-toggle buttons (header + sidebar)
  document.addEventListener("click", function (e) {
    var btn = e.target.closest("button.theme-toggle");
    if (!btn) return;
    e.preventDefault();
    e.stopImmediatePropagation();
    toggle();
  }, true);

  // On load: resolve "auto" to explicit light/dark
  var stored = localStorage.getItem("theme");
  if (!stored || stored === "auto") {
    applyTheme(resolveAuto());
  }
})();
