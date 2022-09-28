context("output class")

set.seed(2022)
a <- raster::raster(matrix(rnorm(81), 9, 9))

test_that("output is a raster layer", {
  expect_equal(class(strucDiv(a, 3, fun = contrast, 
                              na.handling = na.omit, 
                              rank = FALSE))[1], "RasterLayer")
})


test_that("testing output", {
  
  set.seed(2022)
  n <- 5
  a <- raster::raster(matrix(sample(1:4, size=n^2, replace=TRUE), n, n))
  b1 <- raster::as.matrix(strucDiv(a, 3, fun = contrast, 
                                   na.handling = na.pass, 
                                   rank = FALSE) )
#  dput(b1)  # objects are created via dput()
  
  expect_equal(raster::as.matrix(strucDiv(a, 3, fun = contrast)),
               matrix(c(NA, NA, NA, NA, NA, NA, 1.8, 1.15, 1.2, NA, NA, 2, 
                           1.55, 1.75, NA, NA, 1.45, 1.55, 2.4, NA, NA, NA, NA, NA, NA),
                      5L, 5L) )
  expect_equal(raster::as.matrix(strucDiv(a, 3, fun = contrast, na.handling = na.omit)),
               matrix(c(NA, NA, NA, NA, NA, NA, 1.8, 1.15, 1.2, NA, NA, 2, 
                           1.55, 1.75, NA, NA, 1.45, 1.55, 2.4, NA, NA, NA, NA, NA, NA), 
                      5L, 5L) )
  
  set.seed(2024)
  n <- 5
  mata <- matrix(sample(1:4, size=n^2, replace=TRUE), n, n)
  mata[3,3] <- NA
  a <- raster::raster(mata)
  b1 <- raster::as.matrix(strucDiv(a, 3, fun = contrast, 
                                #   na.handling = na.omit, 
                                   rank = FALSE) )
#  dput(b1)  # objects are created via dput()
  
  expect_equal(raster::as.matrix(strucDiv(a, 3, fun = contrast)),
               matrix(NA_real_,n,n)  )
  expect_equal(raster::as.matrix(strucDiv(a, 3, fun = contrast, na.handling = na.omit)),
               matrix(c(0.15, 1.5, 2.7, 2.35, 1, 0.7, 2.1, 3.2, 3.35, 1.95, 
                           0.9, 1.65, 2, 2.25, 1.5, 1.3, 1.95, 1.4, 1.8, 1.15, 0.95, 1.55, 
                           1.2, 1.2, 0.6), 5L, 5L) )
  
  })

test_that("argument testing", {
  expect_error(strucDiv(a, -3, fun = contrast, na.handling = na.omit, 
                        rank = FALSE)) # ncol not positive
  
  expect_error(strucDiv(a, -3, fun = contrast, na.handling = sum, 
                        rank = FALSE)) # na.handling wrong function
  
  expect_error(strucDiv(a, 3, fun = contrast, na.handling = 'na.omit', 
                        rank = FALSE)) # na.handling not passed as function
  
  
  
})  


