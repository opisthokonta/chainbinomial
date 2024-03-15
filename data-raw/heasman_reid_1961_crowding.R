## code to prepare heasman_reid_1961_crowded data set.


heasman_reid_1961_crowding <- read.table(file = fs::path('data-raw', 'heasman_reid_1961_crowding.txt'), header=TRUE)

usethis::use_data(heasman_reid_1961_crowding, overwrite = TRUE)
