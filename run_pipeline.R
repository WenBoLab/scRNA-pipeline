source("scripts/bootstrap.R")
arguments <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
config_path <- if (length(arguments) && !startsWith(arguments[1], "--")) arguments[1] else "config/config.R"
run_pipeline(config_path, resume = !"--no-resume" %in% arguments)
