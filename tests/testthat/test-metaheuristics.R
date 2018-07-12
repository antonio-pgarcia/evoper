library(evoper)

context("Metaheuristic verification")

test_that("PSO metaheuristic results", {
  skip_on_cran()

  level<- elog.level()
  elog.level("ERROR")

  set.seed(2718282)
  f<- PlainFunction$new(f0.cigar4)
  f$Parameter(name="x1",min=-100,max=100)
  f$Parameter(name="x2",min=-100,max=100)
  f$Parameter(name="x3",min=-100,max=100)
  f$Parameter(name="x4",min=-100,max=100)
  v<- extremize("pso", f)
  results<- v$getBest()

  expect_equal(round(results$x1,3), 0)
  expect_equal(round(results$x2,3), 0)
  expect_equal(round(results$x3,3), 0)
  expect_equal(round(results$x4,3), 0)

  elog.level(level)
})


test_that("SAA metaheuristic results", {
  skip_on_cran()

  level<- elog.level()
  elog.level("ERROR")

  set.seed(2718282)
  f<- PlainFunction$new(f0.cigar4)
  f$Parameter(name="x1",min=-100,max=100)
  f$Parameter(name="x2",min=-100,max=100)
  f$Parameter(name="x3",min=-100,max=100)
  f$Parameter(name="x4",min=-100,max=100)
  v<- extremize("saa", f)
  results<- v$getBest()

  expect_equal(round(results$x1,3), 0)
  expect_equal(round(results$x2,3), 0)
  expect_equal(round(results$x3,3), 0)
  expect_equal(round(results$x4,3), 0)

  elog.level(level)
})


test_that("ACOR metaheuristic results", {
  skip_on_cran()

  level<- elog.level()
  elog.level("ERROR")

  set.seed(1133)
  f<- PlainFunction$new(f0.cigar4)
  f$Parameter(name="x1",min=-100,max=100)
  f$Parameter(name="x2",min=-100,max=100)
  f$Parameter(name="x3",min=-100,max=100)
  f$Parameter(name="x4",min=-100,max=100)
  v<- extremize("acor", f)
  results<- v$getBest()

  expect_equal(round(results$x1,3), 0)
  expect_equal(round(results$x2,3), 0)
  expect_equal(round(results$x3,3), 0)
  expect_equal(round(results$x4,3), 0)

  elog.level(level)
})


test_that("EES1 metaheuristic results", {
  skip_on_cran()

  level<- elog.level()
  elog.level("ERROR")

  set.seed(2718282)
  f<- PlainFunction$new(f0.cigar4)
  f$Parameter(name="x1",min=-100,max=100)
  f$Parameter(name="x2",min=-100,max=100)
  f$Parameter(name="x3",min=-100,max=100)
  f$Parameter(name="x4",min=-100,max=100)
  v<- extremize("ees1", f)
  results<- v$getBest()

  expect_equal(round(results$x1,3), 0)
  expect_equal(round(results$x2,3), 0)
  expect_equal(round(results$x3,3), 0)
  expect_equal(round(results$x4,3), 0)

  elog.level(level)
})


test_that("TABU Search metaheuristic results", {
  skip_on_cran()

  level<- elog.level()
  elog.level("ERROR")

  set.seed(2718282)
  f<- PlainFunction$new(f0.cigar4)
  f$Parameter(name="x1",min=-100,max=100)
  f$Parameter(name="x2",min=-100,max=100)
  f$Parameter(name="x3",min=-100,max=100)
  f$Parameter(name="x4",min=-100,max=100)
  v<- extremize("tabu", f)
  results<- v$getBest()

  expect_equal(round(results$x1,3), 0)
  expect_equal(round(results$x2,3), 0)
  expect_equal(round(results$x3,3), 0)
  expect_equal(round(results$x4,3), 0)

  elog.level(level)
})


