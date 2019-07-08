context("invls")

library(raster)

data(dtm100m,dtm1m)
ls1 <- invls(dtm100m, raster::extent(dtm1m), 180)

test_that("invls output", {
  expect_is(ls1, "RasterLayer")
})
