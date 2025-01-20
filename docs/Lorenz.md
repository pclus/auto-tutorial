# Lorenz system

Bifurcation diagram for the Lorenz system.

Check: 

(1) Nonlinear Dynamics and Chaos, chapter 9. Strogatz (1995).

These are commands for the Python interface of auto-07p.

## Provide files:
* lorenz.f90 :  Encodes the system equations, the corresponding Jacobian, and a function to save maxima and minima of limit cycles.
* c.lorenz : constants of auto-07p.
* lorenz.py: code in python/auto-07p to run the bifurcation diagram.

## Tutorial:

Initializing and finding the fixed points as a function of the parameter r. Thanks to our functions GETUY_MIN and GETUY_MAX defined in the jrnmm.f90 file, auto returns the min and max values of 'x' in the limit cycles.

```
# Use Euler integrator up to steady state to initialize the diagram
init=run('lorenz',IPS=-2,NMX=100000,PAR={'r':0.0})                 
ic=init(201)                                                
fp=run(ic,IPS=1,NMX=10000,ISW=1,ICP=[1,2],UZSTOP={'r' : 25.0})
# Continuation along r, from 0 to 25, give us: 
# i) For r < 1 the origin is stable.
# ii) at r = 1 the origin loses its stability through a pitchfork bifurcation, branded 'BP', and two symmetric stable branches emerge.
# iii) at r ~ 21 the branches lose their stability in a subcritical Hopf bifurcation 'HB'.                                            
```
Check that we can access the content of fp by typing:
```
print(fp)
  BR    PT  TY  LAB       r          L2-NORM          x             y             z           x_min     
   1    14  BP  202   1.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  -1.00000E+01
   1   254  UZ  203   2.50000E+01   0.00000E+00   0.00000E+00   0.00000E+00   0.00000E+00  -1.00000E+01

  BR    PT  TY  LAB       r          L2-NORM          x             y             z           x_min     
   2   328  HB  204   2.14286E+01   2.23392E+01   6.39196E+00   6.39196E+00   2.04286E+01  -1.00000E+01
   2   380  UZ  205   2.50000E+01   2.59230E+01   6.92820E+00   6.92820E+00   2.40000E+01  -1.00000E+01
```
Notice that auto takes care of continuing along the branch that emerges at the pitchfork bifurcation automatically.
It only gives us one due to the symmetric nature of it.
The two branches are organized according to their 'BR' value, so to work with them we need to specify them:
```
plt.plot(fp[0]['r'],fp[0]['x']) #to plot the first branch.
plt.plot(fp[1]['r'],fp[1]['x']) #to plot the second branch.
plt.show()
```

Now we continue from the Hopf bifurcation 'HB' and set (IPS=2). 
This branch dies on a SNIC bifurcation 'LP' around r ~ 11

```
hb1=fp('HB1')                                                                                   
lc1=run(fp('HB1'),IPS=2,ISP=2,ICP=[1,11,2,3],ISW=1,DS='-',STOP=['BP1'])
```
This part of the solution can be plotted normally:
```
plt.scatter(lc1['r'],lc1['x_min'], c=2 * (pt_vals_2(lc1) < 0), cmap=greys, s=0.8)
```
