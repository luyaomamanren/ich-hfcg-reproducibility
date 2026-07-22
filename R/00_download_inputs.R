source("R/_setup.R")

bundled_python <- "C:/Users/32574/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
python <- Sys.getenv("PYTHON_FOR_DOWNLOADS", unset = if (file.exists(bundled_python)) bundled_python else "python")
helper <- file.path(project_root, "tools", "bootstrap_download.py")
status <- system2(
  python,
  c(helper, "--out", file.path(project_root, "data", "raw"),
    "--series", "GSE24265", "GSE163256", "GSE266873", "--workers", "6", "--chunk-mb", "1")
)
if (status != 0) stop("GEO download failed; inspect network access and rerun.")
save_session_info("00_download")
