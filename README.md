# sMGPT

Alejandro Aviles
(avilescervante@gmail.com)

#

Mathematica routines to compute the RSD power spectrum for biased tracers in Modified Gravity models using the SPT/EFT theory of https://arxiv.org/abs/2012.05077

#

Run the the notebooks 0_All.nb in folders code/ and codeNW/. Then, run the file pk_EFT_IR.nb to obtain the IR-resummed EFT 1-loop power spectrum:

pk[k_, mu_, b1_, b2_, bs2_, b3nl_, alpha0_, alpha2_, alpha4, ctilde , alphashot0_, alphashot2_, PshotP_]

and compute its multipoles. 

b1, b2, bs2 and b3nl are the bias parameters. alpha0, alpha2, alpha4 the leading order EFT params, and ctilde the Next-to-leading order EFT param. The total shot noise is given by Pshot = PshotP (alphashot0 + alphashot2 k^2 mu^2). Hence, PshotP is intended to be fixed to the Poissonian noise 1/n, and alphashot0 a correction.   


The input power spectra should be located in input/ where you can find a code PNW.nb to compute the non-wiggle power spectrum.

In the modules code/0_params.m and codeNW/0_params.m you should specify the route to the power spectra with and without wiggles. There you can specify the kind of perturbative kernels used. 



Please contact me for any doubt. 

#

The code is based mainly in

Alejandro Aviles, Georgios Valogiannis, Mario A.Rodriguez-Meza, Jorge L. Cervantes-Cota, Baojiu Li & Rachel Bean, Redshift space power spectrum beyond Einstein-de Sitter kernels, JCAP04(2021)039, https://arxiv.org/abs/2012.05077

It uses SPT kernels obtained from a mapping of the LPT kernels of:

Alejandro Aviles & Jorge L. Cervantes-Cota, Lagrangian perturbation theory for modified gravity, Phys. Rev. D 96, 123526 (2017). https://arxiv.org/abs/2012.05077

You can use and modify this code as long as you cite at least the first paper above.






