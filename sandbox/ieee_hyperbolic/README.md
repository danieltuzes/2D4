# 2D4\::sandbox::ieee_hyperbolic
This project helps to identify the effect of using *non-IEEE compatible hyperbolic function calculation* instead of *IEEE compatible hyperbolic function calculation*. The latter is the default but the first supposed to be faster.

This improvements also eliminates checking for the case of close dislocations for the far fields. It is obvious and known apriori, that distance( [-0.5:0.5), ± i ) cannot be < 10^{-6}, where i is integer.

Measurements on the precision say that
* In the stress field, the media of the difference is 1e-16, the average is 1e-14, and the maximum is 1e-11.
* In the derivated field, the media is 1e-11, the average is 9e-11, and the maximum is 8e-8.

It seems that the non-IEEE compatible hyperbolic function calculation is **precise enough**.

The [speedtest is available](https://github.com/danieltuzes/2D4/wiki/Speedtests#ieee_hyperbolic) in Hungarian and says that one can gain a speed improvement by **180% to 200%**. The smalle5r systems can gain a bit more with this improvement.
