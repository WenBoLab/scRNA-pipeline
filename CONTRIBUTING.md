# Contributing

Use R for analysis and tests. Keep README files and guides in English ASCII.
Preserve raw counts, input provenance and explicit candidate-label status.
Run `Rscript tests/run_tests.R` before proposing analytical changes. For changes
that affect demo results, rerun the demo, verification and R Markdown lesson;
update the validation record and examples with the new R results.

Do not commit credentials, private patient data, installed runtime libraries,
FASTQs, BAMs or the generated full `results/` directory. Synthetic tests must
be identified as synthetic. State which biological conclusions are supported
by independent donors and which results are exploratory cell-level markers.
