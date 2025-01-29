---
layout: default
title: Installation
nav_order: 2
---

{% include mathjax.html %}

# Installing `auto-07p`

In order to install `auto-07p` in Linux from the official repository you can use:

```
mkdir auto-07p
git clone https://github.com/auto-07p/auto-07p auto-07p
cd auto-07p
./configure
make
make install
```

Depending on your system, there might be some conflicts.
I suggest installing a minimal version without the provided plotting tools,
as they require some dependencies that are outdated or conflict with current packages.
To do so, diable them in the configuration step of the previous instructions:

```
./configure --enable-plaut04=no --enable-plaut=no --enable-plaut-qt=no
```

Some other conflicts might appear anyhow, please see the official documentation.

If everything goes according to plan, typing `auto` in a terminal should start the interface.


