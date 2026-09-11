source("scripts/bootstrap.R")
check_packages(c(required_packages(), "testthat", "withr"))
test_results <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
test_table <- as.data.frame(test_results)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)
write_csv(test_table[, setdiff(names(test_table), "result"), drop = FALSE], "results/validation/test_results.csv")
write_json(list(implementation = "R", tests = nrow(test_table), failed_expectations = sum(test_table$failed),
  errored_tests = sum(test_table$error), r_version = as.character(getRversion()),
  scope = "Unit checks plus explicitly synthetic end-to-end, Harmony and donor-level DESeq2 tests"),
  "results/validation/test_summary.json")
