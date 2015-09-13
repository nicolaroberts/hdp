if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  test_that("Package style conforms to linters", {
    lintr::expect_lint_free(linters=lintr::with_defaults(
      trailing_whitespace_linter=NULL, infix_spaces_linter=NULL))
  })
}
