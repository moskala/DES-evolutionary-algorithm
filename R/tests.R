source("des.R")
source("functions/ackley.R")
source("functions/shubert.R")
source("functions/shekel.R")
source("functions/rastr.R")
source("functions/griewank.R")
source("functions/perm0db.R")
source("functions/rothyp.R")
source("functions/zakharov.R")
library(data.table)


test <- function(func, expected_value, n_times, test_name, filename){
  errors <- c()
  time_values <- c()
  for(val in 1:n_times) {
    start.time <- Sys.time()
    res <- func();
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    value <- res$value;
    errors <- append(errors, abs(value - expected_value));
    time_values <- append(time_values, as.numeric(time.taken));
  }
  best_error <- min(errors);
  worst_error <- max(errors);
  mean_error <- mean(errors);
  median_error <- median(errors);
  std_error <- sd(errors);
  
  best_time <- min(time_values);
  worst_time <- max(time_values);
  mean_time <- mean(time_values);
  median_time <- median(time_values);
  std_time <- sd(time_values);
  
  line_error = c(test_name, "error", best_error, worst_error, mean_error, median_error, std_error);
  fwrite(as.data.table(t(as.matrix(line_error))), file = filename, append=TRUE, sep = ";")
  line_time = c(test_name, "time", best_time, worst_time, mean_time, median_time, std_time);
  fwrite(as.data.table(t(as.matrix(line_time))), file = filename, append=TRUE, sep = ";")
}

# Ackley Function min = 0; x* = 0
des_ackeley <- function() {return(des(sample(-30:30, N, replace=TRUE), ackley, lower=-30, upper=30))}

# Shubert Function min = -186.7309
des_schubert <- function() {return(des(sample(-10:10, 2, replace=TRUE), shubert, lower=-10, upper=10))}

# Shekel Function min = -10.5364
des_shekel <- function() {return(des(sample(0:10, 4, replace=TRUE), shekel, lower=0, upper=10))}

# Rastrigin Function min = 0
des_rastrigin <- function() {return(des(sample(-5:5, N, replace=TRUE), rastr, lower=-5.12, upper=5.12))}

# Griewank Function min = 0; x* = 0
des_griewank <- function() {return(des(sample(-600:600, N, replace=TRUE), griewank, lower=-600, upper=600))}

# Bowl-Shaped
# Perm Function 0, d, β min = 0
des_perm <- function() {return(des(sample(-600:600, N, replace=TRUE), perm0db, lower=-600, upper=600))}

# Rotated Hyper-Ellipsoid Function min = 0
des_rotated <- function() {return(des(sample(-65:65, N, replace=TRUE), rothyp, lower=-65.536, upper=65.536))}

# Plate-Shaped
# Zakharov Function min = 0
des_zakharov <- function() {return(des(sample(-5:10, N, replace=TRUE), zakharov, lower=-65.536, upper=65.536))}

