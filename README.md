## Example IRLBA ##

Use [my fork](https://github.com/zdk123/irlba) of the irlba package to call the code from C++/Armadillo.

* [Original source code](https://github.com/bwlewis/irlba)
* [Relevant issue](https://github.com/bwlewis/irlba/issues/37)


## Usage ##

```r
library(devtools)
install_github('zdk123/irlba')
install_github('zdk/ExampleIRLBA')

## generate some rank 1 data
x <- rnorm(50) %*% t(rnorm(50))

## check that outputs are identical
y <- rnorm(50)
head(r_orthog(x, y))
head(irlba:::orthog(y,x))

r_irlba(x, 4)$d
irlba::irlba(x, 4)$d
```
