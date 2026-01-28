# Extracted from test-metaheuristics.R:137

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "evoper", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(evoper)
context("Metaheuristic verification")

# test -------------------------------------------------------------------------
skip_on_cran()
level<- elog.level()
elog.level("ERROR")
set.seed(2718282)
f<- PlainFunction$new(f0.cigar4)
f$Parameter(name="x1",min=-100,max=100)
f$Parameter(name="x2",min=-100,max=100)
f$Parameter(name="x3",min=-100,max=100)
f$Parameter(name="x4",min=-100,max=100)
v<- extremize("ga", f)
