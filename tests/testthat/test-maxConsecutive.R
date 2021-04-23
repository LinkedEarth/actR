test_that("maxConsecutive works numeric", {
  mc <- maxConsecutive(c(1,2,7,5,7,7,7,34),val = 7)
  expect_equal(mc$max, 3)
  expect_equal(mc$index,5:7)
})

test_that("maxConsecutive works for TRUE by default", {
  mc <- maxConsecutive(c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE))
  expect_equal(mc$max, 4)
  expect_equal(mc$index,c(3:6,8:11))
})

test_that("maxConsecutive indexing works above minimum values", {
  mc <- maxConsecutive(c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),gte = 3)
  expect_equal(mc$index,c(8:11))
  mc <- maxConsecutive(c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),gte = 2)
  expect_equal(mc$index,c(5,6,8:11))
  mc <- maxConsecutive(c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,TRUE),gte = 5)
  expect_equal(mc$index,NA)
})
