local({
  project <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  Sys.setenv(
    RENV_PATHS_ROOT = file.path(project, ".renv-cache"),
    RENV_PATHS_CACHE = file.path(project, ".renv-cache", "cache"),
    RENV_PATHS_SANDBOX = file.path(project, ".renv-cache", "sandbox")
  )
  source("renv/activate.R")
})
