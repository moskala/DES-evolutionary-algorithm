source("tests.R")

N <- 20
test_count <- 10
test_filename <- "tests_R_results_20_dim.txt"
fwrite(t(c("name", "expected_value", "res_value", "error", "time")), col.names = FALSE, file = test_filename, sep = ";")
seed_value <- 1234
set.seed(seed_value)

test_ackeley <- test(des_ackeley, 0, test_count, "ackeley", test_filename)
test_schubert <- test(des_schubert, -186.7309, test_count, "schubert", test_filename)
test_shekel <- test(des_shekel, -10.5364, test_count, "shekel", test_filename)
test_rastrigin <- test(des_rastrigin, 0, test_count, "rastrigin", test_filename)
test_griewank <- test(des_griewank, 0, test_count, "griewank", test_filename)
test_perm  <- test(des_perm, 0, test_count, "perm", test_filename)
test_rotated <- test(des_rotated, 0, test_count, "rotated", test_filename)
test_zakharov <- test(des_zakharov, 0, test_count, "zakharov", test_filename)
