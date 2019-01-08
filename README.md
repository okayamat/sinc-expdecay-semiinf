# sinc-expdecay-semiinf
Sinc approximation for exponentially decaying functions over the semi-infinite interval

## Overview
These programs approximate the following three functions:

(1) f(t) = t^(pi/4) exp(-t)  
(2) f(t) = (exp(t) - 1)^(1/2) exp(- 3t/2)  
(3) f(t) = (1 + (1 - 2 exp(-t))^2)^(1/2) exp(-t) t / (1 + t)

Approximation is done by means of the Sinc approximation combined with
a proper variable transformation. Programs that begin with Stenger use
a variable transformation given by Stenger [1], t = arcsinh(exp(x)),
whereas other programs use another transformation t = log(1 + exp(x)),
which is described in [2].

Each program investigates maximum approximation error among selected
sampling points increasing n as n = 2, 7, 12, 17, 22, ..., then outputs
n, maximum error, and its error bound. The error bound for Stenger's
approximation is shown in [3], and the error bound for another approximation
is shown in [2].

## Results
Outputs by those programs are stored in data/ directory, with .dat extension.
Gnuplot programs and created eps graphs are also stored in the directory.

computation environment:

OS: Mac OS X 10.9  
CPU: 1.7 GHz Intel Core i7  
Memory: 8 GB 1600 MHz DDR3  
Compiler: Apple LLVM version 6.0  

## References
1. F. Stenger:
 Numerical Methods Based on Sinc and Analytic Functions, Springer-Verlag,
 New York, 1993.
2. T. Okayama:
 New conformal map for the Sinc approximation for exponentially decaying
 functions over the semi-infinite interval, arXiv:1812.11546 [math.NA],
 December 2018.
3. T. Okayama:
 Error estimates with explicit constants for the Sinc approximation over
 infinite intervals, Applied Mathematics and Computation, Vol. 319 (2018),
 pp. 125--137.
