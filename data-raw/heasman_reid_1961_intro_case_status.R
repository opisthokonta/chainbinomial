## code to prepare heasman_reid_1961_intro_case_status data set.


heasman_reid_1961_intro_case_status <- read.table(file = fs::path('data-raw', 'heasman_reid_1961_intro_case_status.txt'), header=TRUE)
colnames(heasman_reid_1961_intro_case_status) <- tolower(colnames(heasman_reid_1961_intro_case_status))

usethis::use_data(heasman_reid_1961_intro_case_status, overwrite = TRUE)
