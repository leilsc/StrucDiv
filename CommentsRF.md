# Remarks RF

- make sure the comments of StrucDiv are incorporated here.

- rename StrucDiv2 to StrucDiv (including .Rproj file)

- why imports of tidyr and parallel?

- nonstandard file README.Rmd: maybe an .Rbuildignore?

- there are quite a few notes and warnings when checking the package


## man

- Diversity:  
  *  `\usage` : why different presentation?
  *   to the StrucDiv.  ->  to StrucDiv.  or to the function StrucDiv().
  * `\example`: not meaningful 
   
- StrucDiv/StrucDivNet:

   * `\value`: A raster layer of the same dimension as the input raster layer.  
      or maybe even a bit more explicit.
   * some text indicating the differences of both functions. 
   * Please set to TRUE only for bigger problems. ???
   * arguments: display_progress 
   
- ndvi/ndvi.15gl:  
   * is the `\describe`-construct adequate here??   
   * possibly `\seealso{ndvi.gl15}` and vice versa
   * in `\example` give the code from ndvi to ndvi.gl15
   
- why is there a function sub-directory?   



## data:

- i would change the content of `ndvi@file@name` and `ndvi.15gl@file@name`


## R


Strucdiv:

- I believe the `return(out)` is too early (if `filename` is given)


StrucDivNest:

- there is a spurious `return()` on line 615.

