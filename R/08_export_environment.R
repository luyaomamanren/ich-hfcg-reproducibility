source("R/_setup.R")
ip <- as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")])
ip <- ip[order(ip$Package), ]
write_tsv(ip, file.path(project_root, "results", "tables", "R_package_versions.tsv"))
capture.output(sessionInfo(), file = file.path(project_root, "results", "logs", "sessionInfo_complete.txt"))
writeLines(R.version.string, file.path(project_root, "results", "logs", "R_version.txt"))
