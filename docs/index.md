---
layout: default
title: Index
nav_order: 1
nav_exclude: false
---

{% include mathjax.html %}

# Introduction

This is a tutorial for numerical bifurcation analysis using the well established **auto-07p** software.
From the many capabilities of this tool, in this tutorial we focus on the analysis of nonlinear dynamical systems.
The tutorial consists of a complete guide-through of several examples of varying levels of complexity. 

The official [**auto-07p**](https://github.com/auto-07p/auto-07p) package also includes a documentation with several demonstrations. 
This tutorial aims to complement those examples.

## What is auto-07p?

AUTO is a numerical continuation software originally developped by Eusebius Doedel in 1976.
The software has been updated in several new realises, from which auto-07p, released in 2007, is the last official version.

The core version of auto-07p is written, as its predecessors, in Fortran. However, auto-07p includes a Python interface which largely simplifies its usage. 

Auto-07p allows for bifurcation analysis of systems of coupled ODEs with initial conditions, boundary value problems, and even parabollic PDEs. Nonetheless, here we focus on the first of these cases.

Therefore, given a continuous dynamical system, auto-07p allows us to explore:

- Contiunation of fixed points and their stability as one system parameter varies.
- Continuation of periodic orbit and their stability as one parameter varies.
- Detection of bifurcations of fixed points (saddle-node, pitchfork, Hopf, ...) and periodic orbits (period doubling, saddle-node of limit-cycles, torus bifurcation, ...) as one parameter of the system is changed.
- Continuation of these bifurcations in a 2-parameter space.
- Continuation of periodic orbits with a fixed period in 2-parameter space.

<!--
## Why auto-07p?

--!>
