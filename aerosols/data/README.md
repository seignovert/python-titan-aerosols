Optical index database
=======================
![License MIT](https://img.shields.io/badge/license-MIT-green.svg)

Tholins Doose
---------------
Filterring of single scattering albedo from [Doose et al (2016)](https://dx.doi.org/10.1016/j.icarus.2015.09.039) complemented with CIRS measurement

> __Note:__
> The values are extrapolated as constant for wavelengths
> above 800 &micro;m and below 20 nm.
>
> In between, the real part is linear interpolated and
> the imaginary part is LOG interpolated.

Tholins CVD
------------
Laboratory produced optical index constant for tholins aerosols based on [Tomasko et al, 2008](https://dx.doi.org/10.1016/j.pss.2007.11.019).

> __Note:__
> The values are extrapolated as constant for wavelengths
> above 314 &micro;m and below 20 nm.
>
> In between, the real part is linear interpolated and
> the imaginary part is LOG interpolated.
>
> Between 935 nm and 1.5 um, the bump of the imaginary part is
> removed and fixed at 7.19e-3.

Authors:
---------
- Pascal Rannou ([Univ. Reims](https://planeto.univ-reims.fr/rannou/)), original data
- Benoît Seignovert ([Univ. Nantes](https://planeto.univ-reims.fr/seignovert/)), database convertion
