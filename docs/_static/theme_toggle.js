document.addEventListener("DOMContentLoaded", function () {
  var btn = document.querySelector("button.theme-toggle");
  if (!btn) return;
  // Replace Furo's 3-state toggle with 2-state (light ↔ dark)
  btn.addEventListener("click", function (e) {
    e.stopImmediatePropagation();
    var current = document.body.dataset.theme;
    var next = (current === "dark") ? "light" : "dark";
    document.body.dataset.theme = next;
    localStorage.setItem("theme", next);
  }, true);
  // Initialize: resolve "auto" to explicit light/dark
  if (!localStorage.getItem("theme") || localStorage.getItem("theme") === "auto") {
    var isDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    var initial = isDark ? "dark" : "light";
    document.body.dataset.theme = initial;
    localStorage.setItem("theme", initial);
  }
});
