# 2D4::sandbox::ieee_hyperbolic
This project helps to identify the effect of using *non-IEEE compatible hyperbolic function calculation* instead of *IEEE compatible hyperbolic function calculation*. The latter is the default but the first supposed to be faster.

Measurements say that
* In the stress field, the media of the difference is 1e-16, the average is 1e-14, and the maximum is 1e-11.
* In the derivated field, the media is 1e-11, the average is 9e-11, and the maximum is 8e-8.

It seems that the non-IEEE compatible hyperbolic function calculation is **precise enough**.