
load(file = "data/examples.RData")

examples <- purrr::map2(
  examples,
  names(examples),
  ~ parse_dCVnet_input(data = .x$x, y = .x$y, family = .y)
)


examples_merged <- purrr::map2(
  examples,
  names(examples),
  ~ dCVnet:::impy_dat_merger(x = .x$x_mat, y = .x$y, family = .y)
)

examples_unmerged <- purrr::map(examples_merged,
                                dCVnet:::impy_dat_unmerger)

test_that("y-merging and unmerging is lossless", {
  expect_equal(purrr::map(examples_unmerged, "x"),
               purrr::map(examples, "x_mat"),
               ignore_attr = TRUE)
  expect_equal(purrr::map(examples_unmerged, "y"),
               purrr::map(examples, "y"),
               ignore_attr = TRUE)
})
