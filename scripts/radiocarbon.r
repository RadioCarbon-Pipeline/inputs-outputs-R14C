library(rcarbon)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
no_trailing_args <- commandArgs(trailingOnly = FALSE)
script_name_index <- which(grepl("^--file=", no_trailing_args))
script_name <- sub("^--file=", "", no_trailing_args[script_name_index])
script_name <- sub("\\.r$", "", basename(script_name))

print(script_name)

available_curves <- c('intcal20', 'intcal13', 'intcal13nhpine16', 'shcal20', 'shcal13', 'shcal13shkauri16', 'marine13', 'marine20')

confidence_interval <- 0.95
step <- 5
time_left <- 10000
time_right <- 0
curve <- "intcal20"
value <- "Everything"

get_value <- function(i) {
  v <- strsplit(config[[i, 1]], "=")[[1]][2]
  return(v)
}

if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).csv", call. = FALSE)
} else if (length(args) == 1) {
  stop("CAREFUL: No output file specified.", call = FALSE)
} else if (length(args) == 2) {
  print("No config file specified. Using default values. (Confidence Interval -> 0.95, Step -> 5 years, No Filtering)")
} else if (length(args) == 3) {
  print("All arguments specified.")
  config <- read.delim2(args[[3]], header = FALSE, sep = "\n")
  step_cfg <- get_value(1)
  confidence_interval_cfg <- get_value(2)
  time_left_cfg <- get_value(3)
  time_right_cfg <- get_value(4)
  curve_cfg <- get_value(5)

  if (is.na(step_cfg)) {
    print("No Step Interval specified, using default value")
  } else {
    step <- as.numeric(step_cfg)
  }

  if (is.na(confidence_interval_cfg)) {
    print("No Confidence Interval specified, using default value")
  } else {
    confidence_interval <- as.numeric(confidence_interval_cfg)
  }

  if (is.na(time_left_cfg)) {
    print("No start time for the time interval specified, using default value")
  } else {
    time_left <- as.numeric(time_left_cfg)
  }

  if (is.na(time_right_cfg)) {
    print("No end time for the time interval specified, using default value")
  } else {
    time_right <- as.numeric(time_right_cfg)
  }

  if (is.na(curve_cfg)) {
    print("No curve specified, using default value")
  } else {
    curve <- curve_cfg
  }

  if (!(curve %in% available_curves)) {
    stop("ERROR: Curve provided is not available! It should be one of: 'intcal20', 'intcal13', 'intcal13nhpine16', 'shcal20', 'shcal13', 'shcal13shkauri16', 'marine20' or 'marine13'.")
  }

  if (time_left <= time_right) {
    stop("ERROR: The left bound for the time interval should be bigger than the right bound!")
  }

  print(paste("Confidence Interval ->", confidence_interval, ", Step ->", step, "years"))
  print(paste("Time interval -> [", time_left, ", ", time_right, "]"))
  print(paste("Curve -> ", curve))

  print("Starting...")
}

file_name <- sub("\\.csv$", "", basename(args[1]))

date <- format(Sys.time(), "%d-%m-%Y@%H:%M:%S")

# determining the separator of the CSV file
first_line <- readLines(args[[1]], n = 1)
comma_count <- length(gregexpr(",", first_line)[[1]])
semicolon_count <- length(gregexpr(";", first_line)[[1]])

if (comma_count > semicolon_count) {
  sep <- ","
} else {
  sep <- ";"
}

c <- read.csv(args[[1]], sep = sep, stringsAsFactors = FALSE)

if (length(args) == 3) {
  column <- get_value(6)
  value <- get_value(7)

  if (is.na(column) || is.na(value)) {
    print("No subsetting applied")
  } else {
    print(paste("Filtering on column:", column, "with value:", value))

    c <- subset(c, c[[column]] == value)
  }
}

if (nrow(c) == 0) {
  stop("CAREFUL: No values match the subsetting provided.")
}

original_col_len <- ncol(c)

c.caldates <- calibrate(x = c$C14Age, errors = c$C14SD, calCurves = "intcal20", eps = 1e-5, ncores = 4, type = "full")

DK.spd <- spd(c.caldates, timeRange = c(time_left, time_right))

pdf(paste(file_name, script_name, "spd", date, "pdf", sep = "."))

plot(DK.spd)
plot(DK.spd, runm = 200, add = TRUE, type = "simple", col = "darkorange", lwd = 1.5, lty = 2) # using a rolling average of 200 years for smoothing

dev.off()

# only consider values that are inside the confidence interval
considered <- function(lst) {
  len <- length(lst[[1]])
  low <- floor(len * (1 - confidence_interval))
  high <- ceiling(len * confidence_interval)

  return(lst$calBP[low:high])
}

reduce_lists <- function(lst) {
  mins <- sapply(lst, function(x) min(considered(x)))
  maxs <- sapply(lst, function(x) max(considered(x)))

  min <- min(mins)
  max <- max(maxs)

  return(list(min, max))
}

mm <- reduce_lists(c.caldates$grids)

len <- length(c.caldates$grids) # number of rows
num_steps <- ((mm[[2]] - mm[[1]]) / step) + 1 # number of new columns

step_values <- seq(from = mm[[1]], to = mm[[2]], by = step)

new_cols <- matrix(0, nrow = len, ncol = num_steps)
colnames(new_cols) <- paste0(step_values)

col <- 1
for (i in seq(from = mm[[1]], to = mm[[2]], by = step)) {
  for (j in 1:len) {
    elems <- c.caldates$grids[[j]]

    indices <- (which(elems$calBP >= i - step / 2 & elems$calBP <= i + step / 2))
    if (length(indices) > 0) {
      tmp <- list()
      for (k in 1:length(indices)) {
        tmp <- c(elems$PrDens[[indices[[k]]]], tmp)
      }
      new_cols[j, col] <- sum(unlist(tmp))
    }
  }
  col <- col + 1
}

mean_values <- numeric(len)  
median_values <- numeric(len)  

for (j in 1:len) {
  row_values <- new_cols[j, ]
  total_weight <- sum(row_values)

  weighted_mean <- sum(row_values * step_values) / sum(row_values)
  mean_values[j] <- weighted_mean

  cumulative_weights <- cumsum(row_values)
  median_index <- which(cumulative_weights >= total_weight / 2)[1] # Find the index where cumulative weight exceeds or equals half the total weight
  weighted_median <- step_values[median_index]
  median_values[j] <- weighted_median
}

CalibratedDates <- c(paste("Probabilities of", c$C14ID))

c$WeightedMean <- mean_values
c$Median <- median_values
c <- cbind(c, CalibratedDates)
c <- cbind(c, new_cols)

extra_columns = 3

cumulative_values <- apply(new_cols, 2, sum)

new_row <- c(rep("", original_col_len + extra_columns), cumulative_values)

new_row[original_col_len + extra_columns] <- "Cumulative Probabilities"

new_row_df <- as.data.frame(t(new_row))
colnames(new_row_df) <- colnames(c)

c <- rbind(c, new_row_df)

# Output csv
write.csv(c, args[[2]], row.names = FALSE)
print(paste(args[[2]], "created."))

s <- sum(cumulative_values)
weight_cumulative_values <- cumulative_values / s

# Violin pdf
pdf(paste(file_name, script_name, "violin", date, "pdf", sep = "."))

# Calculate the density
kde <- density(step_values, weights = weight_cumulative_values)
kde_df <- data.frame(
  Density = kde$y,
  Years = kde$x
)

# Create the violin plot
ggplot(kde_df, aes(x = "", y = Years, weight = Density)) +
  geom_violin(scale = "area", fill = "lightgreen") +
  scale_y_continuous(limits = c(0, 10000)) +
  labs(x = "Density", y = "Years") +
  geom_boxplot(width = 0.05, fill = "white") +
  xlab(value) +
  ylab("Calendar years BP") +
  ggtitle(label = "Distribution of C14 per site", subtitle = "Illustrates periods of biomass burning") +
  theme_bw()
dev.off()

