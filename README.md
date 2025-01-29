# A tutorial for auto-07p

## Technical notes:

In order to establish the dollar sign $ as a LaTeX delimiter in MathJax, we specified the option at `_includes/mathjax.html`, which is then called by all the markdown files.

Therefore, all markdown files should start with the header

```
---
layout: default
title: Title to appear on the right-hand-side menu
nav_order: 3 (a number specifying the order in the menu)
---

{% include mathjax.html %}

# This is a header...


```
