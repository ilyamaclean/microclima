context("coastalTps")

library(raster)

data(tas,landsearatios)
temp <- tas[,,1]
sst <- 10.665
dT <- if_raster(sst - temp, dtm1km)

lsw <- landsearatios[,,7]  # upwind
lsa <- apply(landsearatios, c(1, 2), mean) # mean, all directions
lsw <- if_raster(lsw, dtm100m)
lsa <- if_raster(lsa, dtm100m)

dTf <- coastalTps(dT, lsw, lsa)

test_that("coastalTps output", {
  expect_is(dTf, "RasterLayer")
  expect_equal(raster::res(lsw), raster::res(dTf))
  expect_equal(raster::extent(lsw), raster::extent(dTf))
})
