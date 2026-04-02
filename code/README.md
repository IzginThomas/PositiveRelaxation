# Numerical Experiments

This directory contains code to reproduce the numerical experiments described
in the manuscript.

This code is developed with Matlab version R2025b. To reproduce the
results, start Matlab in this directory and execute the following commands in
the Matlab command window to create the figures and tables shown in the paper.

The reference solutions are produced using 

```bash
ODESolverTestSuite();
```
where the respective test case is set in line 10. E.g., for the Porous Medium Equation, write

```bash
test ='PME'; % choose from lotka_volterra, lin_ad, IsothEuler (only implemented for MPRK22), strat, PME
```
In addition, the relaxation algorithm is turned on and off in line 13:
```bash
relaxation_flag = 1; % 1="on", 0="off"
```
In any case, choose the respective solver for the relaxation algorithm in line 14, by, e.g.
```bash
relax_meth = 'newton'; % newton, regula, bisection, secant. For Dissipative problems always use 'newton'
```

If you wish to turn on adaptivity, remove the first three dots in line 15, which currently reads

```bash
input = {...'adap',{'opt',{1e-3, 1e-3},'PID'},... 
```
Finally, the methods are selected in the lines 25-27, e.g.line 26 reads
```bash
solvers = {{'MPRK43I', [input(:)', 'parameter', [1/2, 3/4]]}};
```



