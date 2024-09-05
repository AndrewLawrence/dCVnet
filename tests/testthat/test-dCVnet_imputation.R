
load(file = "data/examples.RData")

examples <- mapply(
  function(.x, .y) parse_dCVnet_input(data = .x$x, y = .x$y, family = .y),
  examples,
  names(examples),
  SIMPLIFY = FALSE
)

examples_merged <- mapply(function(.x, .y) {
  impy_dat_merger(x = .x$x_mat,
                  y = .x$y,
                  family = .y)
},
examples,
names(examples),
SIMPLIFY = FALSE)

examples_unmerged <- lapply(examples_merged,
                            impy_dat_unmerger)

test_that("y-merging and unmerging is lossless", {
  expect_equal(lapply(examples_unmerged, `[[`, "x"),
               lapply(examples, `[[`, "x_mat"),
               ignore_attr = TRUE)
  expect_equal(lapply(examples_unmerged, `[[`, "y"),
               lapply(examples, `[[`, "y"),
               ignore_attr = TRUE)
})
