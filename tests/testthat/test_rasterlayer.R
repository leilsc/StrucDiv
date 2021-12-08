context("output class")


a <- raster::raster(matrix(rnorm(81), 9, 9))

test_that("output is a raster layer", {
  expect_equal(class(StrucDiv(a, 3, fun = contrast, 
                              na.handling = na.omit, 
                              rank = FALSE))[1], "RasterLayer")
})


