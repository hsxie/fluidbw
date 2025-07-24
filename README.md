# fluidbw
Linear Fluid Plasma Model with Accurate Kinetic Bernstein Waves

Although kinetic models provide high physical fidelity, they are computationally intensive due to the requirement of evaluating velocity-space integrals. In contrast, fluid models dramatically reduce computational cost by taking velocity moments of the kinetic equation, yielding multi-fluid and magnetohydrodynamic (MHD) formulations. However, these conventional fluid models inherently fail to capture key kinetic effects such as Landau damping and Bernstein modes.

This study demonstrates the first successful construction of a fluid model that captures both kinetic Landau damping and Bernstein modes through precise rational approximations and transformation into an equivalent matrix representation. The fluid matrix method avoids numerical difficulties such as divergence and the need for good initial values when solving the kinetic dispersion relation. This approach can be directly applied to plasma wave propagation conditions analysis and ray tracing studies, and we anticipate further developments in fast full-wave simulations of ICRF using fluid models with accurate kinetic effects. 

Ref: H. S. Xie, Developing a Linear Fluid Plasma Model with Accurate Kinetic Bernstein Waves: A First Step, 2025, https://arxiv.org/abs/2502.06550

-em3d_Lpole.m, The fluid matrix method code to solve kperp. 

-bonk/.., The conventional root finding codes to solve kperp.


Huasheng XIE (huashengxie@gmail.com, ENN)
10:20 2025/7/24
