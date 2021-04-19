test_that("maxConsecutive works numeric", {
  mc <- maxConsecutive(c(1,2,7,5,7,7,7,34),val = 7)
  expect_equal(mc$max, 3)
  expect_equal(mc$index,5:7)
})

test_that("maxConsecutive works for TRUE by default", {
  mc <- maxConsecutive(c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE))
  expect_equal(mc$max, 4)
  expect_equal(mc$index,3:6)
})
