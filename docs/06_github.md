# Upload to flee70973-coder/F

Extract `F-github-upload.zip`. Upload the extracted project contents to the
repository root, where `README.md` should be visible immediately. The ZIP is a
transport package; uploading only the ZIP does not create a browsable project.

Include `README.md`, `LICENSE`, `run_pipeline.R`, the RStudio project,
`renv.lock`, `environment.yml`, `VALIDATION.md`, the manifest, and all project
folders (`R`, `config`, `scripts`, `notebooks`, `docs`, `resources`, `examples`,
`tests`, `renv`, and `.github`). Include `.Rprofile`, `.gitattributes` and
`.gitignore`. Do not upload runtime libraries or newly generated `results/`.

The archive contains the required small teaching data and example results.
Large RDS intermediate objects are delivered separately in the full-results
archive and are not required for the repository.

If replacing the earlier version, remove old Python files, `src/`, `Snakefile`,
`pyproject.toml`, `uv.lock`, the old `.ipynb` lesson and obsolete Python CI files.
Overlaying an R archive does not automatically delete old remote files. The new
archive is a complete replacement repository. Preserve any unrelated personal
files before removing obsolete project content.

Use GitHub's **Add file > Upload files**, drag the extracted contents and commit.
If your operating system hides dotfiles, enable their display so CI and renv
activation are included. GitHub Desktop is another option for replacing files
and committing a directory without writing Git commands.

After uploading, open **Actions**. `R checks` runs the R validation suite.
`PBMC3k R demo` is a manual workflow that runs the real-data demo, verifies raw
counts and resume behavior, renders the lesson and uploads its outputs as a
workflow artifact. These remote workflows can only be assessed after GitHub
successfully receives the project.
