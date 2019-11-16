# 2D4\::sandbox::ieee_hyperbolic
This project helps to identify the effect of using *non-IEEE compatible hyperbolic function calculation* instead of *IEEE compatible hyperbolic function calculation*. The latter was implemented first but implementation of the first supposed to be faster. It turns out that this is indeed the case and the precision is acceptable.

This improvements also eliminates checking for the case of close dislocations for the far fields. It is obvious and known apriori, that distance( [-0.5:0.5), ± i ) cannot be < 10^{-6}, where i is integer.

Further improvements with regard to this line is implemented not investigated here, namely that instead of calculating `sinh` and `cosh` one can calculate only `exp` and derive the required function from it.

## Hyperbolic identities
The first non-IEEE compatible way to calculate cosh and sinh functions is applied on the calculation of cosh[2π * (x±i)] and on sinh[2π * (x±i)] by using the [identities of the hyperbolic function of their argument sum](https://en.wikipedia.org/wiki/Hyperbolic_function#Sums_of_arguments). This makes it possible to calculate only cosh(2π) and sinh(2π), and all the other calculations, cosh[2π * (x±i)] and sinh[2π * (x±i)], i∈{1, 2, 3, 4}, can be applied. E.g.:
```
cosh[2π(x±1)] = cosh(2πx) cosh(2π) ± sinh(2πx) sinh(2π)
sinh[2π(x±1)] = sinh(2πx) cosh(2π) ± cosh(2πx) sinh(2π)
```

Measurements on the precision say that
* In the stress field, the median of the difference is 1e-16, the average is 1e-14, and the maximum is 1e-11.
* In the derivated field, the median is 1e-11, the average is 9e-11, and the maximum is 8e-8.

## Hyperbolic definitions and identities
Realising that not only cosh(x) but sinh(x) is also required for the calculation, and that cosh(x) is a more expensive calculation, one can arrive to the desire to derive the value of cosh(x) from sinh(x), or even further, to derive the value of sinh(x) and cosh(x) from the same quantity. Luckily, it is the case, as both of them can be expressed with the exponential function exp(x) as:
```
cosh(x) = [exp(x) + exp(-x)] / 2
sinh(x) = [exp(x) - exp(-x)] / 2
```
Furthermore, using the hyperbolic identities, for **i⪰3** one can write that sinh(2πi) = cosh(2πi) and factoring this term simplifies the expression to
```
cosh[2π(x±i)] =   cosh(2πi) exp(±2πx)
sinh[2π(x±i)] = ± cosh(2πi) exp(±2πx)
```
Therefore, only exp(x) is needed to be calculated, exp(-x) = 1/exp(x) can be derived with a division, and the hyperbolic functions can be constructed as an addition or substraction and a division.

## Measurements on the precision
By calculating the absolute value of the relative difference between the IEEE-compatible calculation and the speed-improved ones, one can get the precision of these methods. Although on the functions themselves, the error is smaller then on the field and on its derivative field, the latter will be used, therefore, the precision calculation is mode on them too. Below can be found 2 table showing the summary of the statistics of the two non-IEEE compatible calculations.

| only identity | median |  mean | maximum |
|:-------------:|:------:|:-----:|:-------:|
|  stress field |  1e-16 | 1e-14 |  1e-11  |
| derived field |  1e-11 | 9e-11 |   8e-8  |

| defs and identity | median |  mean | maximum |
|:-----------------:|:------:|:-----:|:-------:|
|    stress field   |  2e-16 | 3e-14 |  3e-10  |
|   derived field   |  3e-9  |  4e-8 |   3e-4  |

Using only the relative difference, and not their absolute value tells that the meadian is half as large and negative, so is the mean for the stress field. For the derived field, it is 1/4 and 1/6, and also negative.

## Conclusion

It seems that both non-IEEE compatible hyperbolic function calculation is **precise enough**. To further analyse the results and to measure the effect of this improvement, I will investigate the relaxed systems for the two case with size N=1024.


The [speedtest is available](https://github.com/danieltuzes/2D4/wiki/Speedtests#ieee_hyperbolic) in Hungarian and says that one can gain a speed improvement by **180% to 200%**. The smalle5r systems can gain a bit more with this improvement.