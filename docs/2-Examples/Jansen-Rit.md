---
layout: default
title: Jansen-Rit
nav_order: 5
---

{% include mathjax.html %}


# Jansen-Rit


The Jansen's, or Jansen-Rit, neural mass model was developed in the 90s to model
the emergence of oscillatory activity in cortical columns. 
A detailed description of the biophysical grounds of the system,
the governing equations, and their dynamics can be found in the seminal paper of 
Grimbert and Faugeras: 

<a id="1">[1]</a> 
F. Grimbert, O. Faugeras; Bifurcation Analysis of Jansen's Neural Mass Model. Neural Comput 2006; 18 (12): 3052–3068 [(link)](https://doi.org/10.1162/neco.2006.18.12.3052)

Here we provide the necessary files and a steb-by-steb commands for the Python interface of auto-07p
to generate the bifurcation diagram of the system.

## Provided files:

* `jrnmm.f90`:  Encodes the system equations, the corresponding Jacobian, and a function to save maxima and minima of limit-cycles in the system variable $y=y_1-y_2$. [(download it here)](https://github.com/pclus/auto-tutorial/blob/main/src/Jansen/jrnmm.f90)
* `c.jrnmm`: constants of auto-07p. [(download it here)](https://github.com/pclus/auto-tutorial/blob/main/src/Jansen/c.jrnmm)

## Tutorial:

We start by obtaining the fixed points as a function of the parameter $p$.
Thanks to the functions `GETUY_MIN` and `GETUY_MAX` defined in the `jrnmm.f90` file, auto-07p returns the values of 'y' as defined in [1].

```
# Use Euler integrator up to steady state to initialize the diagram
init=run('jrnmm',IPS=-2,NMX=100000,PAR={'p':-100.0}) # Start from p = -100.
ic=init(201)
# Continue along p, up to p = 500, to uncover the bifurcations (2 saddle-node and 3 Hopf)                                                
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=['p'],UZSTOP={'p' : 500.0})
```

Notice that we are using `ISW = 1`,  which means that we are carrying a one-dimensional bifurcation exploration. Also, we set our continuation parameter to be $p$ with
 `ICP=['p']` (this would be equivalent to set `ICP=[1]`).
Check that we can access the content of `fp` by typing:
```
print(fp)
 BR    PT  TY  LAB       p          L2-NORM          y0            y1            y2            y3            y4            y5      
   1   500      202  -5.03669E+01   2.72358E+00   7.73097E-04  -1.01358E+00   2.52795E+00   0.00000E+00   0.00000E+00   0.00000E+00
   1  1000      203  -3.95059E-01   2.66490E+00   1.90701E-03   6.64140E-01   2.58081E+00   0.00000E+00   0.00000E+00   0.00000E+00
   1  1500      204   4.95732E+01   3.65034E+00   4.69670E-03   2.43955E+00   2.71543E+00   4.84676E-29   0.00000E+00   0.00000E+00
   1  2000      205   9.95211E+01   5.60723E+00   1.23470E-02   4.65895E+00   3.12010E+00   1.59943E-27   0.00000E+00   0.00000E+00
   1  2143  LP  206   1.13586E+02   7.20558E+00   2.08702E-02   6.21929E+00   3.63874E+00  -1.69414E-23   0.00000E+00   0.00000E+00
   1  2500      207   7.80770E+01   9.67214E+00   3.54597E-02   8.44185E+00   4.72061E+00   1.94140E-28   0.00000E+00   0.00000E+00
   1  3000      208   2.81048E+01   1.13359E+01   4.49815E-02   9.86651E+00   5.58141E+00   1.80072E-23   0.00000E+00   0.00000E+00
   1  3500      209  -2.18546E+01   1.33204E+01   5.54463E-02   1.15172E+01   6.69205E+00  -4.75156E-28   0.00000E+00   0.00000E+00
   1  3698  LP  210  -4.13014E+01   1.55772E+01   6.61077E-02   1.33511E+01   8.02452E+00  -4.23326E-24   0.00000E+00   0.00000E+00
   1  3994  HB  211  -1.21475E+01   1.89401E+01   7.98955E-02   1.60292E+01   1.00888E+01   4.71226E-27   0.00000E+00   0.00000E+00
...
```

Auto has detected four bifurcations: two saddle-node (`LP`) and three Hopf (`HB`).
The fixed points can be easily plotted with `matplotlib.pyplot`:
```
import matplotlib.pyplot as plt
plt.plot(fp['p'],fp['y1']-fp['y2'])
plt.show()
```
The stability of each fixed point is stored in the `PT` column, where negative values represent stability, and positive values indicate instability. 
We can easily access them by defining the function:
```
import numpy as np
def pt_vals(f):
	return np.array([f[0][i]['PT'] for i in range(len(f[0]))])

```
We can use this function to plot the fixed points with their stability:
```
from matplotlib import colormaps
plt.scatter(fp['p'],fp['y1'] - fp['y2'], c=2 * (pt_vals(fp) < 0), cmap='Purples', s=0.1)
```
To display the bifurcation points in our plot, we can define the function:
```
def bifs(f,par):
	exceptions = ['No Label', 'RG', 'EP', 'UZ']  # List of exceptions
	return [[f[0][i]['TY name'],f[0][i][par]] for i in range(len(f[0])) if f[0][i]['TY name'] not in exceptions]

```
and add it to our plot:

```
bfp = bifs(fp,'p')
for b in bfp:
	plt.axvline(b[1],color="black", ls="--",alpha=0.7)

plt.scatter(fp['p'],fp['y1'] - fp['y2'], c=2 * (pt_vals(fp) < 0), cmap='Purples', s=0.1)
```


Now we can proceed to compute the limit-cycles and their stability.

We start from the first Hopf bufrcation and set `IPS=2`.
This branch dies on a SNIC bifurcation, characterized by an infinite period and thus we stop the continuation when the period reaches 100.0.
Auto always stores the period of a limit-cycle in the auxiliary parameter 11, which is named `'PERIOD'`:

```
hb1=fp('HB1')                                                                                   
lc1=run(hb1,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, UZSTOP={'PERIOD' : 100.0}) 
```

The second Hopf produces a limit cycle that vanishes at another Hopf (which auto labels as `BP` when continuing limit-cycles)

```
hb2=fp('HB2')
lc2=run(hb2,IPS=2,ISP=2,ICP=[1,11,2,3,4],NMX=20000,ISW=1, STOP=['BP1'])
```



