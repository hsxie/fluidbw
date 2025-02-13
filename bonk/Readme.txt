This is fsolve version to solve the EM3D Maxwellian kinetic dispersion relation (KDR).

- main_tst.m, main program to run.  Input other parameters to solve kperp (kx). You can easily modify it to solve omega with input kz&kx.

-fDR.m, Ronnmark's method to fast calculate KDR, which is not suitable for less interesting strong damped mode, but fast and convergent for all sum_n. Use J-pole for zfun, and RYLA.m to fast calculate Bessel sum_n.

- fD.m & fDtensor.m, conventional method to calculate KDR, with truncate sum_n to n_max, which also use Gamman.m, Gammanp.m and zfun.m to calculate the Gamma_n and Z functions.

Huasheng XIE (huashengxie@gmail.com, ENN)
8:20 2025/2/13