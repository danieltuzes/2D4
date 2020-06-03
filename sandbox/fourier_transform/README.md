# 2D4\::sandbox::fourier_transform
Simple tutorial how to make 2D Fourier transformation

## Description
This program calcualtes the 2D Fourier transformation of the given data. It can do forward and backward Fourier transformation too, based on FFTW. This is not an efficient program but a demonstration tool.

## File format
The data must be in the following format:
```
z(x_0,y_0)	z(x_1,y_0)	z(x_2,y_0)	..........
z(x_0,y_1)	z(x_1,y_1)	z(x_2,y_1)	..........
z(x_0,y_2)	z(x_1,y_2)	z(x_2,y_2)	..........
..........	..........	..........	..........
```

Output data will be presented in both the same format as input format and in 
```
x_0	y_0	z(x_0,y_0)
x_1	y_0	z(x_0,y_0)
...	...	..........

x_0	y_1	z(x_0,y_0)
x_1	y_1	z(x_0,y_0)
...	...	..........

...	...	..........
...	...	..........
...	...	..........
```
### r2c
For real to complex, z **input** value is real. If the size of date is n_x × n_y, then the size of the **output** will be n_x/2 + 1 × n_y.

### c2r
For complex to real, the **input** value z is complex, and the real and imaginary part should be separated with a whitespace (no comma or parenthesis). Include only the nontrivial part of the 2D data. To get a real result, the data must be Hermitian, that is,
$$z(x_i,y_j) = z(x_{n_x-i},y_j)^*$$

where n_x is the size of the data in x direction and \* is the complex conjugation. Therefore, if n_x values in the x is available, feed the program only n_x/2+1 of them.

The size of output will be n_x × n_y, but don't forget to feed only n_x/2 + 1 element in the x direction.