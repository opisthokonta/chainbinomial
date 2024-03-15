## code to prepare heasman_reid_1961_chains data set.



heasman_reid_1961_chains <- read.table(file = fs::path('data-raw', 'heasman_reid_1961_chains.txt'))
colnames(heasman_reid_1961_chains) <- c('chain', 'n')


usethis::use_data(heasman_reid_1961_chains, overwrite = TRUE)
