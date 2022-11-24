test_that("right number of points at vertices", {
  expect_equal(npoints(makepos(spatstat.data::simplenet,0)),
               npoints(vertices.linnet(spatstat.data::simplenet)))
})

test_that("right number of points at vertices in case with duplicates", {
  expect_equal(npoints(makepos(spatstat.data::simplenet,0,duplicate=TRUE)),
               2*nsegments(spatstat.data::simplenet))
})
