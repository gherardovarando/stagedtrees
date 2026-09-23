## data where the outcome does not depend on the treatment at all, so an
## unconstrained search has every reason to merge the two arms
arm_data <- function(n = 2000, seed = 5) {
  set.seed(seed)
  X <- sample(c("x1", "x2"), n, TRUE)
  TT <- sample(c("a", "b"), n, TRUE)
  p <- ifelse(X == "x1", 0.2, 0.8)          # Y depends on X only
  Y <- ifelse(runif(n) < p, "yes", "no")
  data.frame(X = factor(X), TT = factor(TT), Y = factor(Y))
}

test_that("best_move_cpp refuses a move into a stage holding the same group", {
  ct <- matrix(c(10, 90, 12, 88, 50, 50, 55, 45), nrow = 4, byrow = TRUE)
  asg <- c(0L, 1L, 2L, 3L)
  ## situations 0 and 1 are the two arms of one context, 2 and 3 of another
  ctx <- c(0L, 0L, 1L, 1L)
  free <- best_move_cpp(ct, asg, 4L, 0)
  held <- best_move_cpp(ct, asg, 4L, 0, ctx)
  ## unconstrained, the obvious move is to put the two similar ones together
  expect_true(any(free[, 1] %in% c(1, 2) & free[, 2] %in% c(1, 2)))
  ## constrained, no move may join two situations of one context
  joins <- held[held[, 2] > 0, , drop = FALSE]
  for (r in seq_len(nrow(joins))) {
    s <- joins[r, 1]
    target <- joins[r, 2]
    expect_false(ctx[s] %in% ctx[asg == (target - 1)])
  }
})

test_that("stages_hc with separate keeps the arms apart", {
  d <- arm_data()
  f <- full(d, lambda = 1)
  expect_equal(nrow(arm_separation(f, "TT", "Y")), 0)  # full() is separated

  free <- stages_hc(f)
  held <- stages_hc(f, separate = "TT")
  ## the unconstrained search merges the arms, since Y ignores TT
  expect_gt(nrow(arm_separation(free, "TT", "Y")), 0)
  expect_equal(nrow(arm_separation(held, "TT", "Y")), 0)
})

test_that("a constrained search beats repairing an unconstrained one", {
  d <- arm_data()
  f <- full(d, lambda = 1)
  held <- stages_hc(f, separate = "TT")
  repaired <- separate_arms(stages_hc(f), "TT", "Y")
  expect_equal(nrow(arm_separation(repaired, "TT", "Y")), 0)
  ## both keep the arms apart, but the constrained search chose where to
  ## spend the stages instead of having them split after the fact
  expect_lte(BIC(held), BIC(repaired))
})

test_that("separate only constrains the variables which follow the treatment", {
  d <- arm_data()
  f <- full(d, lambda = 1)
  ## X precedes TT, so its staging is searched as usual
  expect_equal(
    stages_hc(f, scope = "X", separate = "TT")$stages$X,
    stages_hc(f, scope = "X")$stages$X
  )
})

test_that("stages_hc refuses a starting model whose arms are already together", {
  d <- arm_data()
  m <- stages_hc(full(d, lambda = 1))
  expect_gt(nrow(arm_separation(m, "TT", "Y")), 0)
  expect_error(stages_hc(m, separate = "TT"), "apart")
  ## and accepts it once repaired
  expect_silent(stages_hc(separate_arms(m, "TT", "Y"), separate = "TT"))
})

test_that("separate leaves the unconstrained search untouched", {
  d <- arm_data()
  f <- full(d, lambda = 1)
  expect_equal(stages_hc(f, separate = NULL)$stages, stages_hc(f)$stages)
})
