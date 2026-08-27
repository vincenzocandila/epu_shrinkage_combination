

############################################################
# Replication notes
############################################################

# The state-level stock market data used in the paper are subject
# to copyright restrictions and therefore cannot be distributed
# with the replication package.
#
# The complete empirical analysis requires the file:
#
#   close_new.csv
#
# Users without access to these proprietary data can nevertheless
# reproduce all results based exclusively on publicly available data.
# In particular, the following parts of the script can be reproduced:
#
#   - Table A.1: state-level EPU summary statistics;
#   - Table A.2: neighborhood-based EPU statistics
#                (columns 1-14; RV statistics are excluded);
#   - Table A.3: summary statistics for U.S. EPU,
#                Industrial Production, and Housing Starts.
#
# Table 1, the RV columns of Table A.2, model estimation,
# forecast combinations, and Tables 3-7 require access to
# close_new.csv.
#
# The script automatically checks whether close_new.csv is available.
# Sections requiring the proprietary stock market data should be
# executed only when this file is present.


##################### let us load the libraries

library(rumidas)		# Version 0.1.3
library(roll)			# Version 1.2.1
library(maxLik)		# Version 1.5.2.2
library(xts)			# Version 0.14.2
library(rugarch)		# Version 1.5.6
library(highfrequency)# Version 1.0.3
library(xtable)		# Version 1.8.8
library(np)			# Version 0.70.5
library(datasets)		# Version 4.6.1
library(lubridate)	# Version 1.9.5
library(fBasics)		# Version 4052.98
library(DescTools) 	# Version 0.99.60
library(readxl)		# Version 1.5.0
library(splm)			# Version 1.6.5
library(caret)		# Version 7.0.1
library(glmnet)		# Version 5.0
library(DT)			# Version 0.34.0

############################################################
# U.S. state spatial adjacency matrix
############################################################

# Load the U.S. state adjacency matrix included in the splm package.
# The object "usaww" provides the spatial weights used to identify
# neighboring U.S. states in the empirical analysis.

data("usaww", package = "splm")

adj_mat <- usaww

##################### let us load the functions 

source(file="comb_functions.R")

##################################################
################################################## MODEL spec
################################################## 

spec_gjr_n <- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)), 
		mean.model=list(armaOrder=c(0,0), include.mean=FALSE),  
		distribution.model="norm")


############################################################
# Order U.S. states by resident population
############################################################

# States are ordered according to their total resident population.
# Population data are obtained from the U.S. Census Bureau:
# https://www2.census.gov/programs-surveys/popest/tables/2020-2023/state/detail/SCPRC-EST2023-18+POP.xlsx
#
# The Excel file must be downloaded from the link above and placed
# in the working directory before running this script.

if (!file.exists("SCPRC-EST2023-18+POP.xlsx")) {
  stop(
    paste0(
      "'SCPRC-EST2023-18+POP.xlsx' was not found.\n",
      "Please download the U.S. Census Bureau population file from:\n",
      "https://www2.census.gov/programs-surveys/popest/tables/2020-2023/",
      "state/detail/SCPRC-EST2023-18+POP.xlsx\n",
      "and place it in the working directory."
    )
  )
}

db_pop <- as.data.frame(
  read_excel(
    "SCPRC-EST2023-18+POP.xlsx",
    skip = 8
  )
)

# Remove the leading "." from state names.
db_pop[, 1] <- sub("^\\.", "", db_pop[, 1])

# Keep only the state name and total resident population.
# The first 51 rows contain the 50 states plus the District of Columbia.
db_pop <- db_pop[1:51, 1:2]

colnames(db_pop) <- c(
  "State",
  "Total Resident Population"
)

# Order states from the largest to the smallest population.
ordered_states <- db_pop[
  order(db_pop[, 2], decreasing = TRUE),
]

# Exclude the District of Columbia, since the analysis considers
# the 50 U.S. states only.
ordered_states <- ordered_states[
  ordered_states[, 1] != "District of Columbia",
  1,
  drop = FALSE
]

ordered_states <- as.character(unlist(ordered_states))

############################################################
# Reorder the spatial adjacency matrix
############################################################

# Convert state names to the format used in the row and column
# names of the adjacency matrix provided by the splm package.
ordered_states_up <- toupper(ordered_states)
ordered_states_up <- gsub(" ", "_", ordered_states_up)

# In the usaww object, Tennessee is labeled "TENNESSE".
ordered_states_up[ordered_states_up == "TENNESSEE"] <- "TENNESSE"

rownames(adj_mat)

############################################################
# Import state-level stock market indices
############################################################

# The file "close_new.csv" contains the daily closing prices of the
# 50 U.S. state-level stock market indices used in the paper.
#
# These data are subject to copyright restrictions and therefore
# cannot be redistributed with the replication package.
#
# Users who have access to the original data should place the file
# "close_new.csv" in the working directory before running this script.
#
# The expected structure is:
#   - first column: Dates (format dd/mm/yyyy);
#   - remaining columns: daily closing prices for the 50 U.S. states,
#     ordered alphabetically by state name

if (!file.exists("close_new.csv")) {
  stop(
    paste0(
      "'close_new.csv' was not found.\n",
      "This file contains copyrighted state-level stock market data ",
      "and cannot be distributed with the replication package.\n",
      "Please obtain the original data from the data provider and ",
      "place 'close_new.csv' in the working directory."
    )
  )
}

db_close <- read.csv(
  file = "close_new.csv",
  header = TRUE,
  dec = ",",
  sep = ";"
)

dim(db_close)

############################################################
# Prepare state-level stock market data
############################################################

# Identify and remove duplicated observations.
# A row is considered duplicated when the closing prices of all
# 50 state-level indices are identical to those of a previous row.
db_close$duplicated <- ifelse(
  duplicated(db_close[, 2:ncol(db_close)]),
  1,
  0
)

db_close <- subset(
  db_close,
  db_close$duplicated == 0
)

# Remove the temporary duplicate indicator.
db_close <- db_close[, -ncol(db_close)]

dim(db_close)


############################################################
# Convert data to an xts object
############################################################

# Convert the first column to dates.
Date_db <- strptime(
  db_close[, 1],
  "%d/%m/%Y",
  tz = "GMT"
)

# Convert state-level closing prices to an xts object.
db_close_i <- as.xts(
  db_close[, 2:ncol(db_close)],
  Date_db
)

head(db_close_i)
tail(db_close_i)
dim(db_close_i)
colnames(db_close_i)


############################################################
# Assign state names and order states by population
############################################################

# In the original file, the 50 state-level indices are stored
# in alphabetical order. Therefore, the column names can be
# assigned using the built-in vector datasets::state.name.
colnames(db_close_i) <- datasets::state.name

N <- ncol(db_close_i)

head(db_close_i)

# Reorder the 50 states from the largest to the smallest
# resident population.
db_close_i <- db_close_i[, ordered_states]

############################################################
# Import state-level Economic Policy Uncertainty data
############################################################

# State-level Economic Policy Uncertainty (EPU) data are obtained from:
# https://www.policyuncertainty.com/state_epu.html
#
# Direct download of the Excel file:
# https://www.policyuncertainty.com/media/State_Policy_Uncertainty.xlsx
#
# The dataset contains three monthly EPU measures for each U.S. state:
# EPU_National, EPU_State, and EPU_Composite.
#
# Download the Excel file and place it in the working directory
# before running this script.

if (!file.exists("State_Policy_Uncertainty.xlsx")) {
  stop(
    paste0(
      "'State_Policy_Uncertainty.xlsx' was not found.\n",
      "Please download the state-level EPU data from:\n",
      "https://www.policyuncertainty.com/media/",
      "State_Policy_Uncertainty.xlsx\n",
      "and place the file in the working directory."
    )
  )
}

epu <- as.data.frame(
  read_excel("State_Policy_Uncertainty.xlsx")
)

############################################################
# Construct monthly dates
############################################################

epu$month_2 <- ifelse(
  epu$month < 10,
  paste0("0", epu$month),
  epu$month
)

epu$month_2 <- paste0(epu$month_2, "/")
epu$day <- "01/"

Dates <- paste0(
  epu$day,
  epu$month_2,
  epu$year
)

# The empirical sample used in the paper ends in December 2024.
# Since the source file is periodically updated, observations after
# December 2024 are excluded to ensure reproducibility of the results.

keep <- as.Date(Dates, format = "%d/%m/%Y") <= as.Date("2024-12-01")

epu   <- epu[keep, ]
Dates <- Dates[keep]

dim(epu)
head(epu)
tail(epu)

# Remove the original year and month variables, together with
# the temporary variables used to construct the monthly dates.
epu <- epu[
  , !names(epu) %in% c("year", "month", "month_2", "day")
]

############################################################
# Reshape state-level EPU data
############################################################

TT <- nrow(epu)

epu_list <- list()

# The dataset contains three EPU measures for each of the
# 50 U.S. states plus the District of Columbia.
for (i in seq(1, 153, by = 3)) {

  col_sel <- c(i, i + 1, i + 2)
  j <- col_sel[3] / 3

  Name <- colnames(epu[, col_sel])[1]

  # Extract the two-letter state abbreviation from the
  # name of the first EPU variable.
  State_Abb <- substr(
    Name,
    nchar(Name) - 1,
    nchar(Name)
  )

  epu_list[[j]] <- cbind(
    Dates,
    State_Abb,
    epu[, col_sel]
  )

  colnames(epu_list[[j]])[3:5] <- c(
    "EPU_National",
    "EPU_State",
    "EPU_Composite"
  )
}


############################################################
# Remove the District of Columbia
############################################################

epu_list <- epu_list[
  sapply(epu_list, function(x) x$State_Abb[1] != "DC")
]

length(epu_list)
# 50

head(epu_list[[1]])
tail(epu_list[[1]])

# Stack all state-level datasets.
epu <- do.call(rbind, epu_list)

############################################################
# Convert EPU data to state-specific xts objects
############################################################

epu_by_state_i <- list()

for (i in 1:50) {

  x <- subset(
    epu,
    epu$State_Abb == state.abb[i]
  )

  Date_x <- strptime(
    x$Dates,
    "%d/%m/%Y",
    tz = "GMT"
  )

  epu_by_state_i[[i]] <- as.xts(
    x[, 3:5],
    Date_x
  )
}

names(epu_by_state_i) <- state.name

head(epu_by_state_i[[1]])

# Reorder states according to resident population.
epu_by_state_i <- epu_by_state_i[ordered_states]


############################################################
# Import U.S. Economic Policy Uncertainty data
############################################################

# U.S. Economic Policy Uncertainty (EPU) data are obtained from:
# https://www.policyuncertainty.com/us_monthly.html
#
# Direct download of the Excel file:
# https://www.policyuncertainty.com/media/US_Policy_Uncertainty_Data.xlsx
#
# Download the Excel file and place it in the working directory
# before running this script.
#
# Note that the source dataset is periodically updated or revised.
# The empirical results reported in the paper are based on the version
# of the dataset downloaded in March 2025. The version currently
# available from the source differs slightly from the one originally
# used in the paper and may therefore produce small differences in
# the corresponding summary statistics and empirical results.

if (!file.exists("US_Policy_Uncertainty_Data.xlsx")) {
  stop(
    paste0(
      "'US_Policy_Uncertainty_Data.xlsx' was not found.\n",
      "Please download the U.S. EPU data from:\n",
      "https://www.policyuncertainty.com/media/",
      "US_Policy_Uncertainty_Data.xlsx\n",
      "and place the file in the working directory."
    )
  )
}

epu_full <- as.data.frame(
  read_excel("US_Policy_Uncertainty_Data.xlsx",
sheet="Legacy Three Component EPU")
)

epu_full<-epu_full[complete.cases(epu_full),]

############################################################
# Construct monthly dates
############################################################

epu_full$Month_2 <- ifelse(
  epu_full$Month < 10,
  paste0("0", epu_full$Month),
  epu_full$Month
)

epu_full$Month_2 <- paste0(epu_full$Month_2, "/")
epu_full$Day <- "01/"

epu_full$Date <- paste0(
  epu_full$Day,
  epu_full$Month_2,
  epu_full$Year
)

# The empirical sample used in the paper ends in December 2024.
# Since the source file is periodically updated, observations after
# December 2024 are excluded to ensure reproducibility of the results.
keep <- as.Date(
  epu_full$Date,
  format = "%d/%m/%Y"
) <= as.Date("2024-12-01")

epu_full <- epu_full[keep, ]

Date_epu_full <- strptime(
  epu_full$Date,
  "%d/%m/%Y",
  tz = "GMT"
)

# Retain the Three-Component Index, which is the U.S. EPU measure
# used in the empirical analysis.
epu_full_i <- as.xts(
  epu_full[, "Three_Component_Index"],
  Date_epu_full
)

tail(epu_full_i)

################################################## 
################################################## OTHER MIDAS VARIABLES
################################################## (first realise)

############################################################
# Import vintage Industrial Production data
############################################################

# Vintage data for the U.S. Industrial Production Index (INDPRO)
# are obtained from ALFRED:
# https://alfred.stlouisfed.org/series/downloaddata?seid=INDPRO
#
# Units:
# Levels
#
# In the download page, select vintage dates from:
# 1990-01-17 to 2025-01-17
#
# Output format:
# "Observations by Vintage Date, All Observations"
#
# File format:
# Excel
#
# Save the downloaded file as:
# "vintage_industrial_production.xlsx"
#
# Place the file in the working directory
# before running this script.

if (!file.exists("vintage_industrial_production.xlsx")) {
  stop(
    paste0(
      "'vintage_industrial_production.xlsx' was not found.\n",
      "Please download the INDPRO vintage data from ALFRED and ",
      "place the file in the working directory."
    )
  )
}

ip <- as.data.frame(
  read_excel(
    "vintage_industrial_production.xlsx",
    sheet = "Vintages Starting 1990-01-17"
  )
)

############################################################
# Restrict the observation sample
############################################################

# The downloaded ALFRED file contains historical observations
# starting in 1919. The empirical analysis only requires
# observations from January 1990 onward.
ip <- ip[
  as.Date(ip[, 1]) >= as.Date("1990-01-01"),
]

############################################################
# Construct the real-time Industrial Production series
############################################################

ip_i <- data.frame(
  matrix(
    NA,
    ncol = 2,
    nrow = nrow(ip)
  )
)

colnames(ip_i) <- c("Date", "IP")

ip_i[, 1] <- ip[, 1]

# Remove the observation-date column.
ip <- ip[, -1]

NCOL <- ncol(ip)

# For each observation date, retain the earliest vintage in which
# the Industrial Production observation is available. This produces
# the real-time series used in the empirical analysis.
for (i in 1:nrow(ip)) {

  ip_which <- which(
    !is.na(ip[i, 1:NCOL])
  )[1]

  ip_i[i, 2] <- ip[i, ip_which]
}

############################################################
# Convert Industrial Production data to xts
############################################################

Date_ip <- strptime(
  ip_i$Date,
  "%Y-%m-%d",
  tz = "GMT"
)

ip_i <- as.xts(
  ip_i[, 2],
  Date_ip
)

plot(ip_i)

# Compute the annualized monthly growth rate.
ip_f_i <- ((ip_i / lag(ip_i))^12 - 1) * 100

ip_f_i[1] <- 0

############################################################
# Import vintage Housing Starts data
############################################################

# Vintage data for U.S. Housing Starts (HOUST)
# are obtained from ALFRED:
# https://alfred.stlouisfed.org/series/downloaddata?seid=HOUST
#
# Units:
# Thousands of Units
#
# In the download page, select vintage dates from:
# 1990-01-18 to 2025-01-17
#
# Output format:
# "Observations by Vintage Date, All Observations"
#
# File format:
# Excel
#
# Save the downloaded file as:
# "vintage_housing_starts.xlsx"
#
# Place the file in the working directory
# before running this script.

if (!file.exists("vintage_housing_starts.xlsx")) {
  stop(
    paste0(
      "'vintage_housing_starts.xlsx' was not found.\n",
      "Please download the Housing Starts vintage data from ALFRED ",
      "and place the file in the working directory."
    )
  )
}

hs <- as.data.frame(
  read_excel(
    "vintage_housing_starts.xlsx",
    sheet = "Vintages Starting 1990-01-18"
  )
)

############################################################
# Restrict the observation sample
############################################################

# The downloaded ALFRED file contains historical observations
# prior to the period needed for the empirical analysis.
# Retain observations from January 1990 onward.
hs <- hs[
  as.Date(hs[, 1]) >= as.Date("1990-01-01"),
]

############################################################
# Construct the real-time Housing Starts series
############################################################

hs_i <- data.frame(
  matrix(
    NA,
    ncol = 2,
    nrow = nrow(hs)
  )
)

colnames(hs_i) <- c("Date", "HS")

hs_i[, 1] <- hs[, 1]

# Remove the observation-date column.
hs <- hs[, -1]

NCOL <- ncol(hs)

# For each observation date, retain the earliest vintage in which
# the Housing Starts observation is available. This produces the
# real-time series used in the empirical analysis.
for (i in 1:nrow(hs)) {

  hs_which <- which(
    !is.na(hs[i, 1:NCOL])
  )[1]

  hs_i[i, 2] <- hs[i, hs_which]
}

############################################################
# Convert Housing Starts data to xts
############################################################

Date_hs <- strptime(
  hs_i$Date,
  "%Y-%m-%d",
  tz = "GMT"
)

hs_i <- as.xts(
  hs_i[, 2],
  Date_hs
)

plot(hs_i)

# Compute the annualized monthly growth rate.
hs_f_i <- ((hs_i / lag(hs_i))^12 - 1) * 100

hs_f_i[1] <- 0

summary(hs_f_i)

plot(hs_f_i)

############################################################
# Table 1: Summary statistics by state
############################################################

# Table 1 reports summary statistics for state-level stock returns,
# with states ordered from the largest to the smallest resident population.
#
# The starting date differs across states because the corresponding
# state-level stock market indices have different data availability.
#
# Replication of the published values requires access to the copyrighted
# state-level stock market data contained in "close_new.csv".

############################################################
# Define state groups according to data availability
############################################################

# Five state indices use data from January 2010 onward.
states_2010 <- ordered_states[c(33, 40, 45, 46, 50)]

# Seven state indices use data from October 2003 onward.
states_2003 <- ordered_states[c(14, 32, 35, 38, 43, 48, 49)]

# All remaining state indices use data from February 2001 onward.

############################################################
# Initialize Table 1
############################################################

table_1 <- as.data.frame(
  matrix(
    NA,
    ncol = 8,
    nrow = 50
  )
)

############################################################
# Compute state-level summary statistics
############################################################

for (i in 1:50) {

  state_i <- ordered_states[i]

  ##########################################################
  # Select the state-specific sample
  ##########################################################

  if (state_i %in% states_2010) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2010-01-05/2024"]

  } else if (state_i %in% states_2003) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2003-10-02/2024"]

  } else {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2001-02-28/2024"]
  }

  ##########################################################
  # Winsorize and annualize returns
  ##########################################################

  # Winsorize daily returns at the 0.1% and 99.9% quantiles.
  r_t <- Winsorize(
    r_t,
    val = quantile(
      r_t,
      probs = c(0.001, 0.999),
      na.rm = FALSE
    )
  )

  # Express daily returns in annualized percentage units.
  r_t <- r_t * 100 * sqrt(252)

  ##########################################################
  # Compute summary statistics
  ##########################################################

  N_i    <- length(r_t)
  min_i  <- min(r_t)
  max_i  <- max(r_t)
  mean_i <- mean(r_t)
  sd_i   <- sd(r_t)
  skew_i <- as.numeric(skewness(r_t))
  kurt_i <- as.numeric(kurtosis(r_t))

  table_1[i, 1] <- state_i

  table_1[i, 2:8] <- c(
    N_i,
    min_i,
    max_i,
    mean_i,
    sd_i,
    skew_i,
    kurt_i
  )
}

############################################################
# Assign column names
############################################################

colnames(table_1) <- c(
  "State",
  "Obs.",
  "Minimum",
  "Maximum",
  "Mean",
  "Std. Dev.",
  "Skewness",
  "Kurtosis"
)

table_1

############################################################
# Table A.1: Summary statistics of state-level EPUs
############################################################

# Table A.1 reports summary statistics for the three state-level
# EPU measures:
#
#   - EPU_National
#   - EPU_State
#   - EPU_Composite
#
# together with their state-specific average ("Avg.").
#
# All data used to construct this table are publicly available
# from https://www.policyuncertainty.com/state_epu.html
#
# The sample period differs across states because of differences
# in the availability of the corresponding state-level stock
# market indices used in the empirical analysis.
#
# The resulting number of monthly observations is:
#
#   - 300 observations: January 2000 - December 2024
#   - 280 observations: September 2001 - December 2024
#   - 204 observations: January 2008 - December 2024
#
# States are reported in decreasing order of resident population.

############################################################
# Initialize Table A.1
############################################################

table_a.1 <- as.data.frame(
  matrix(
    NA,
    ncol = 18,
    nrow = 50
  )
)

############################################################
# Compute state-level summary statistics
############################################################

for (i in 1:50) {

  state_i <- ordered_states[i]

  ##########################################################
  # Select the state-specific sample
  ##########################################################

  if (state_i %in% states_2010) {

    time_midas <- "2008-01/2024"

  } else if (state_i %in% states_2003) {

    time_midas <- "2001-09/2024"

  } else {

    time_midas <- "2000-01/2024"
  }

  ##########################################################
  # Select state-level EPU data
  ##########################################################

  X_var <- epu_by_state_i[[i]][time_midas]

  ##########################################################
  # National EPU component
  ##########################################################

  X_var_nat <- X_var[, "EPU_National"]

  X_var_nat <- X_var_nat[
    complete.cases(X_var_nat),
  ]

  X_var_nat <- X_var_nat[2:length(X_var_nat)]

  # Replace zero values before computing log differences.
  X_var_nat[
    X_var_nat == 0
  ] <- mean(
    X_var_nat[X_var_nat != 0]
  )

  # Express the monthly log-difference in annualized
  # percentage terms.
  X_var_nat <- makeReturns(X_var_nat) * 100 * 12

  ##########################################################
  # State EPU component
  ##########################################################

  X_var_state <- X_var[, "EPU_State"]

  X_var_state <- X_var_state[
    complete.cases(X_var_state),
  ]

  X_var_state <- X_var_state[2:length(X_var_state)]

  # Replace zero values before computing log differences.
  X_var_state[
    X_var_state == 0
  ] <- mean(
    X_var_state[X_var_state != 0]
  )

  # Express the monthly log-difference in annualized
  # percentage terms.
  X_var_state <- makeReturns(X_var_state) * 100 * 12

  ##########################################################
  # Composite EPU component
  ##########################################################

  X_var_composite <- X_var[, "EPU_Composite"]

  X_var_composite <- X_var_composite[
    complete.cases(X_var_composite),
  ]

  X_var_composite <- X_var_composite[
    2:length(X_var_composite)
  ]

  # Replace zero values before computing log differences.
  X_var_composite[
    X_var_composite == 0
  ] <- mean(
    X_var_composite[X_var_composite != 0]
  )

  # Express the monthly log-difference in annualized
  # percentage terms.
  X_var_composite <- makeReturns(
    X_var_composite
  ) * 100 * 12

  ##########################################################
  # Average of the three EPU components
  ##########################################################

  X_var_all_mean <- apply(
    X_var,
    1,
    mean
  )

  X_var_all_mean <- as.xts(
    X_var_all_mean,
    time(X_var)
  )

  X_var_all_mean <- X_var_all_mean[
    complete.cases(X_var_all_mean),
  ]

  X_var_all_mean <- X_var_all_mean[
    2:length(X_var_all_mean)
  ]

  # Replace zero values before computing log differences.
  X_var_all_mean[
    X_var_all_mean == 0
  ] <- mean(
    X_var_all_mean[X_var_all_mean != 0]
  )

  # Express the monthly log-difference in annualized
  # percentage terms.
  X_var_all_mean <- makeReturns(
    X_var_all_mean
  ) * 100 * 12

  ##########################################################
  # Store summary statistics
  ##########################################################

  table_a.1[i, 1] <- state_i

  # Number of monthly observations in the original EPU sample.
  table_a.1[i, 2] <- nrow(X_var)

  # National EPU: Min, Max, Mean, Std. Dev.
  table_a.1[i, 3:6] <- summ_stat_f(
    X_var_nat
  )

  # State EPU: Min, Max, Mean, Std. Dev.
  table_a.1[i, 7:10] <- summ_stat_f(
    X_var_state
  )

  # Composite EPU: Min, Max, Mean, Std. Dev.
  table_a.1[i, 11:14] <- summ_stat_f(
    X_var_composite
  )

  # Average EPU: Min, Max, Mean, Std. Dev.
  table_a.1[i, 15:18] <- summ_stat_f(
    X_var_all_mean
  )

 }

############################################################
# Assign column names
############################################################

colnames(table_a.1) <- c(
  "State",
  "Obs.",
  "Nat. Min",
  "Nat. Max",
  "Nat. Mean",
  "Nat. Std. Dev.",
  "Loc. Min",
  "Loc. Max",
  "Loc. Mean",
  "Loc. Std. Dev.",
  "Comp. Min",
  "Comp. Max",
  "Comp. Mean",
  "Comp. Std. Dev.",
  "Avg. Min",
  "Avg. Max",
  "Avg. Mean",
  "Avg. Std. Dev."
)

table_a.1

############################################################
# Table A.2: Neighborhood-based EPUs and realized volatility
############################################################

# Table A.2 reports summary statistics for:
#
#   - neighborhood-based National EPU;
#   - neighborhood-based Local (State) EPU;
#   - neighborhood-based Composite EPU;
#   - monthly Realized Volatility (RV).
#
# The three neighborhood-based EPU measures can be fully replicated
# using the publicly available state-level EPU data and the U.S.
# adjacency matrix provided by the splm package.
#
# The RV columns require the copyrighted state-level stock market
# data contained in "close_new.csv" and are therefore constructed
# separately below.
#
# States are reported in decreasing order of resident population.

############################################################
# Initialize Table A.2
############################################################

table_a.2 <- as.data.frame(
  matrix(
    NA,
    ncol = 18,
    nrow = 50
  )
)

############################################################
# Compute neighborhood-based EPU measures
############################################################

for (i in 1:50) {

  state_i <- ordered_states[i]

  ##########################################################
  # Select the state-specific sample
  ##########################################################

  if (state_i %in% states_2010) {

    year_b <- 2008
    period_of_int <- "2008/2024"

  } else if (state_i %in% states_2003) {

    year_b <- 2001
    period_of_int <- "2001/2024"

  } else {

    year_b <- 2000
    period_of_int <- "2000/2024"
  }

  ##########################################################
  # Alaska and Hawaii
  ##########################################################

  # Alaska and Hawaii have no contiguous U.S. neighbors.
  # Therefore, neighborhood-based EPU measures are not defined.
  if (state_i %in% c("Alaska", "Hawaii")) {

    Nat_neigh_mean_f   <- NA
    State_neigh_mean_f <- NA
    Comp_neigh_mean_f  <- NA

    n_obs <- if (state_i %in% states_2010) 204 else
             if (state_i %in% states_2003) 280 else 300

##########################################################
# Maine
##########################################################

} else if (state_i == "Maine") {

  # Maine has only one neighboring U.S. state: New Hampshire.
  # Restrict the neighboring-state EPU series to the sample
  # used for Maine before applying any transformation.
  X_var <- epu_by_state_i[["New Hampshire"]][period_of_int]

  ########################################################
  # National EPU
  ########################################################

  X_var_nat <- X_var[, "EPU_National"]

  X_var_nat <- X_var_nat[
    complete.cases(X_var_nat),
  ]

  X_var_nat <- X_var_nat[
    2:length(X_var_nat)
  ]

  # Replace zero values before computing log differences.
  X_var_nat[X_var_nat == 0] <-
    mean(X_var_nat[X_var_nat != 0])

  Nat_neigh_mean_f <-
    makeReturns(X_var_nat) * 100 * 12

  ########################################################
  # Local EPU
  ########################################################

  X_var_state <- X_var[, "EPU_State"]

  X_var_state <- X_var_state[
    complete.cases(X_var_state),
  ]

  X_var_state <- X_var_state[
    2:length(X_var_state)
  ]

  # Replace zero values before computing log differences.
  X_var_state[X_var_state == 0] <-
    mean(X_var_state[X_var_state != 0])

  State_neigh_mean_f <-
    makeReturns(X_var_state) * 100 * 12

  ########################################################
  # Composite EPU
  ########################################################

  X_var_composite <- X_var[, "EPU_Composite"]

  X_var_composite <- X_var_composite[
    complete.cases(X_var_composite),
  ]

  X_var_composite <- X_var_composite[
    2:length(X_var_composite)
  ]

  # Replace zero values before computing log differences.
  X_var_composite[X_var_composite == 0] <-
    mean(X_var_composite[X_var_composite != 0])

  Comp_neigh_mean_f <-
    makeReturns(X_var_composite) * 100 * 12

  # Maine belongs to the main 2000-2024 sample.
  n_obs <- nrow(X_var)

  ##########################################################
  # All other contiguous states
  ##########################################################

  } else {

    # Identify neighboring states using the adjacency matrix.
    adj_mat_1 <- adj_mat[ordered_states_up[i], ]

    # Add zero entries for Hawaii and Alaska, which are not included
    # in the original usaww adjacency matrix.
    adj_mat_i <- c(adj_mat_1, 0, 0)

    names(adj_mat_i)[49:50] <- c(
      "HAWAII",
      "ALASKA"
    )

    # Reorder according to the population-based state ordering.
    adj_mat_i_f <- adj_mat_i[ordered_states_up]

    # Positions of neighboring states.
    vec_i <- which(adj_mat_i_f > 0)

    list_of_neighborhood <- epu_by_state_i[vec_i]

    ########################################################
    # Construct neighborhood EPU arrays
    ########################################################

    NN <- length(year_b:2024) * 12

    X_var_neigh_array <- array(
      NA,
      dim = c(
        NN,
        length(vec_i),
        3
      )
    )

    for (j in seq_along(vec_i)) {

      X_var_neigh <- list_of_neighborhood[[j]]

      X_var_neigh_array[, j, 1] <-
        X_var_neigh[, "EPU_National"][period_of_int]

      X_var_neigh_array[, j, 2] <-
        X_var_neigh[, "EPU_State"][period_of_int]

      X_var_neigh_array[, j, 3] <-
        X_var_neigh[, "EPU_Composite"][period_of_int]
    }

    X_var_neigh_array[
      is.nan(X_var_neigh_array) |
      is.infinite(X_var_neigh_array)
    ] <- NA

    ########################################################
    # Average EPU measures across neighboring states
    ########################################################

    # Use the dates from one of the state-level EPU series.
    dates_neigh <- time(
      epu_by_state_i[[i]][period_of_int]
    )

    if (length(vec_i) > 1) {

      Nat_neigh_mean <- apply(
        X_var_neigh_array[, , 1],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

      State_neigh_mean <- apply(
        X_var_neigh_array[, , 2],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

      Comp_neigh_mean <- apply(
        X_var_neigh_array[, , 3],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

    } else {

      Nat_neigh_mean <-
        X_var_neigh_array[, , 1]

      State_neigh_mean <-
        X_var_neigh_array[, , 2]

      Comp_neigh_mean <-
        X_var_neigh_array[, , 3]
    }

    Nat_neigh_mean <- as.xts(
      Nat_neigh_mean,
      dates_neigh
    )

    State_neigh_mean <- as.xts(
      State_neigh_mean,
      dates_neigh
    )

    Comp_neigh_mean <- as.xts(
      Comp_neigh_mean,
      dates_neigh
    )

    ########################################################
    # Convert to annualized percentage log-differences
    ########################################################

    Nat_neigh_mean_f <-
      makeReturns(Nat_neigh_mean) * 100 * 12

    State_neigh_mean_f <-
      makeReturns(State_neigh_mean) * 100 * 12

    Comp_neigh_mean_f <-
      makeReturns(Comp_neigh_mean) * 100 * 12

    # Set the first observation to zero, as in the original
    # empirical implementation.
    Nat_neigh_mean_f[1]   <- 0
    State_neigh_mean_f[1] <- 0
    Comp_neigh_mean_f[1]  <- 0

    n_obs <- nrow(Nat_neigh_mean)
  }

  ##########################################################
  # Store publicly replicable results
  ##########################################################

  table_a.2[i, 1] <- state_i
  table_a.2[i, 2] <- n_obs

  if (!state_i %in% c("Alaska", "Hawaii")) {

    table_a.2[i, 3:6] <-
      summ_stat_f(Nat_neigh_mean_f)

    table_a.2[i, 7:10] <-
      summ_stat_f(State_neigh_mean_f)

    table_a.2[i, 11:14] <-
      summ_stat_f(Comp_neigh_mean_f)
  }

  message(i)
}

############################################################
# Assign column names
############################################################

colnames(table_a.2) <- c(
  "State",
  "Obs.",
  "Neigh. Nat. Min",
  "Neigh. Nat. Max",
  "Neigh. Nat. Mean",
  "Neigh. Nat. Std. Dev.",
  "Neigh. Loc. Min",
  "Neigh. Loc. Max",
  "Neigh. Loc. Mean",
  "Neigh. Loc. Std. Dev.",
  "Neigh. Comp. Min",
  "Neigh. Comp. Max",
  "Neigh. Comp. Mean",
  "Neigh. Comp. Std. Dev.",
  "RV Min",
  "RV Max",
  "RV Mean",
  "RV Std. Dev."
)
 
### Table A.2 without the last three columns

table_a.2

############################################################
# Add Realized Volatility statistics to Table A.2
############################################################

# The Realized Volatility (RV) statistics require the copyrighted
# state-level stock market data contained in "close_new.csv".
#
# Therefore, these columns can only be reproduced by users who have
# access to the original state-level stock market data.
#
# The first 14 columns of Table A.2, based on publicly available EPU
# data, can be reproduced independently of this block.

for (i in 1:50) {

  state_i <- ordered_states[i]

  ##########################################################
  # Select the state-specific sample
  ##########################################################

  if (state_i %in% states_2010) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2010-01-20/2024"]

    time_midas <- "2008-01/2024"

  } else if (state_i %in% states_2003) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2003-09-24/2024"]

    time_midas <- "2001-09/2024"

  } else {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2002-01-02/2024"]

    time_midas <- "2000-01/2024"
  }

  ##########################################################
  # Construct monthly Realized Volatility
  ##########################################################

  # Daily returns over the MIDAS sample.
  r_t_rv <- makeReturns(
    db_close_i[, i]
  )[time_midas]

  # Winsorize daily returns using the 0.1% and 99.9% quantiles
  # computed from the corresponding state-specific return sample.
  r_t_rv <- Winsorize(
    r_t_rv,
    val = quantile(
      r_t,
      probs = c(0.001, 0.999),
      na.rm = TRUE
    )
  )

  # Monthly realized volatility, expressed in annualized
  # percentage terms.
  RV_monthly <- sqrt(
    apply.monthly(
      (r_t_rv * 100)^2,
      sum
    ) * 12
  )

  ##########################################################
  # Store RV summary statistics
  ##########################################################

  table_a.2[i, 15:18] <- summ_stat_f(
    RV_monthly
  )

  message(i)
}

table_a.2


############################################################
# Table A.3: Summary statistics
############################################################

# Construct the variables over the common sample 2000-2024.

epu_full_i_subsample <- makeReturns(epu_full_i) * 100 * 12
epu_full_i_subsample <- epu_full_i_subsample["2000/2024"]

ip_f_i_subsample <- ip_f_i["2000/2024"]
hs_f_i_subsample <- hs_f_i["2000/2024"]

# Initialize the table.
table_a.3 <- as.data.frame(
  matrix(
    NA,
    ncol = 6,
    nrow = 3
  )
)

# Compute summary statistics.
table_a.3[1, 3:6] <- summ_stat_f(epu_full_i_subsample)
table_a.3[2, 3:6] <- summ_stat_f(ip_f_i_subsample)
table_a.3[3, 3:6] <- summ_stat_f(hs_f_i_subsample)

# Variable labels and number of observations.
table_a.3[, 1] <- c("Glob.", "IP", "HS")
table_a.3[, 2] <- length(epu_full_i_subsample)

colnames(table_a.3) <- c(
  "",
  "Obs.",
  "Min",
  "Max",
  "Mean",
  "Std. Dev."
)

table_a.3

############################################################
# Cross-state average EPU measures: M5-M7
############################################################

# Construct monthly matrices containing the three EPU measures
# for all 50 states over the common period 1995-2024.

Nat <- State <- Comp <- matrix(
  NA,
  nrow = 360,
  ncol = 50
)

colnames(Nat) <-
  colnames(State) <-
  colnames(Comp) <-
  ordered_states

for (i in 1:50) {

  X_var <- epu_by_state_i[[i]]["1995/2024"]

  Nat[, i]   <- X_var[, "EPU_National"]
  State[, i] <- X_var[, "EPU_State"]
  Comp[, i]  <- X_var[, "EPU_Composite"]
}

############################################################
# M5: Cross-state average National EPU
############################################################

Nat[
  is.nan(Nat) | is.infinite(Nat)
] <- NA

Nat_mean <- apply(
  Nat,
  1,
  function(x) mean(x, na.rm = TRUE)
)

Nat_mean <- as.xts(
  Nat_mean,
  time(X_var)
)

# Express monthly changes in annualized percentage terms.
Nat_mean <- makeReturns(Nat_mean) * 100 * 12


############################################################
# M6: Cross-state average Local EPU
############################################################

State[
  is.nan(State) | is.infinite(State)
] <- NA

State_mean <- apply(
  State,
  1,
  function(x) mean(x, na.rm = TRUE)
)

State_mean <- as.xts(
  State_mean,
  time(X_var)
)

# Express monthly changes in annualized percentage terms.
State_mean <- makeReturns(State_mean) * 100 * 12


############################################################
# M7: Cross-state average Composite EPU
############################################################

Comp[
  is.nan(Comp) | is.infinite(Comp)
] <- NA

Comp_mean <- apply(
  Comp,
  1,
  function(x) mean(x, na.rm = TRUE)
)

Comp_mean <- as.xts(
  Comp_mean,
  time(X_var)
)

# Express monthly changes in annualized percentage terms.
Comp_mean <- makeReturns(Comp_mean) * 100 * 12


############################################################
# Model estimation
############################################################

# All GARCH-MIDAS models are estimated using rumidas::ugmfit().
# The estimation procedure is deterministic and does not rely on
# random-number generation. Therefore, no random seed is required
# to reproduce the estimation results.

# Estimation settings.
#
# K     : number of monthly MIDAS lags
# lstep : forecast horizon and re-estimation frequency, in trading days
# Tin   : rolling in-sample estimation window, in daily observations

K <- 24
lstep <- 25

Tin <- 1000

############################################################
# Initialize result containers
############################################################

# In-sample estimation results for each state.
in_sample_by_state <- list()

# Out-of-sample forecasts for each state.
oos_by_state <- list()

# Multi-step out-of-sample forecasts for each state.
oos_by_state_multi_step <- list()

# Additional container used to store state-specific
# in-sample estimation objects.
in_sample_by_state_list_array <- list()

############################################################
# Initialize table of first-window MIDAS coefficient estimates
############################################################

tab_mat_est <- matrix(
  NA,
  ncol = 15,
  nrow = 50
)

colnames(tab_mat_est) <- c(
  "State-i",
  "Nat.",
  "State",
  "Comp.",
  "Mean",
  "Mean-Nat.",
  "Mean-State",
  "Mean-Comp.",
  "Mean-Nat.-Neigh.",
  "Mean-State-Neigh.",
  "Mean-Comp.-Neigh.",
  "Full",
  "RV",
  "IP",
  "HS"
)

############################################################
# Select states to be estimated
############################################################

# For simplicity, all 50 states are estimated in a single run,
# according to the population-based ordering. Users may split
# the estimation into smaller subsets of states and subsequently
# aggregate the resulting objects.

i_b <- 1
i_e <- 50

############################################################
# Model specifications
############################################################

# m1  : GARCH-MIDAS + state-specific National EPU
# m2  : GARCH-MIDAS + state-specific Local EPU
# m3  : GARCH-MIDAS + state-specific Composite EPU
# m4  : GARCH-MIDAS + average of the three state-specific EPU measures
#
# m5  : GARCH-MIDAS + cross-state average National EPU
# m6  : GARCH-MIDAS + cross-state average Local EPU
# m7  : GARCH-MIDAS + cross-state average Composite EPU
#
# m8  : GARCH-MIDAS + neighboring-state average National EPU
# m9  : GARCH-MIDAS + neighboring-state average Local EPU
# m10 : GARCH-MIDAS + neighboring-state average Composite EPU
#
# m11 : GARCH-MIDAS + U.S. aggregate EPU
# m12 : GARCH-MIDAS + monthly Realized Volatility
# m13 : GARCH-MIDAS + Industrial Production
# m14 : GARCH-MIDAS + Housing Starts
#
# m15 : GJR-GARCH benchmark


############################################################
# State-level estimation loop
############################################################

# Models are estimated separately for each state in the selected
# subset. The full analysis can therefore be split into different
# groups of states and run in separate R sessions.

####################################### 
####################################### BEGIN cycle
#######################################

for (i in i_b:i_e) {

  state_i <- ordered_states[i]

  tab_mat_est[i, 1] <- state_i

  # File containing the final estimation results for this
  # subset of states.
  filename <- paste0(
    "est_new_", i_b, "_", i_e, ".RData"
  )

  # Temporary file used to save intermediate results during
  # the rolling estimation procedure.
  temp_filename <- paste0(
    "temp_new_", i_b, "_", i_e, ".RData"
  )

  message(
    "Estimating state ", i, " of 50: ", state_i
  )

  ############################################################
  # State-level daily returns
  ############################################################

  # The starting date differs across states according to the
  # availability of the corresponding state-level stock index.

  if (state_i %in% states_2010) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2010-01-05/2024"]

  } else if (state_i %in% states_2003) {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2003-10-02/2024"]

  } else {

    r_t <- makeReturns(
      db_close_i[, i]
    )["2001-02-28/2024"]
  }

  # A longer return series is used to construct monthly
  # realized volatility.

  r_t_rv <- makeReturns(
    db_close_i[, i]
  )["1998/2024"]


  ############################################################
  # Winsorize and annualize daily returns
  ############################################################

  # Winsorization thresholds are computed from the state-specific
  # estimation sample and applied both to the estimation returns
  # and to the longer return series used to construct RV.
  winsor_limits <- quantile(
    r_t,
    probs = c(0.001, 0.999),
    na.rm = TRUE
  )

  r_t <- Winsorize(
    r_t,
    val = winsor_limits
  )

  r_t_rv <- Winsorize(
    r_t_rv,
    val = winsor_limits
  )

  # Express daily returns in annualized percentage units.
  r_t <- r_t * 100 * sqrt(252)

  ############################################################
  # Monthly realized volatility
  ############################################################

  # Construct monthly realized volatility from daily returns
  # and express it in annualized percentage terms.
  RV_monthly <- sqrt(
    apply.monthly(
      (r_t_rv * 100)^2,
      sum
    ) * 12
  )


  ############################################################
  # State-specific EPU variables: M1-M4
  ############################################################

  X_var <- epu_by_state_i[[i]]

  # M1: state-specific National EPU
  X_var_nat <- X_var[, "EPU_National"]
  X_var_nat <- X_var_nat[complete.cases(X_var_nat), ]
  X_var_nat <- X_var_nat[2:length(X_var_nat)]
  X_var_nat[X_var_nat == 0] <-
    mean(X_var_nat[X_var_nat != 0])

  X_var_nat <- makeReturns(X_var_nat) * 100 * 12


  # M2: state-specific Local EPU
  X_var_state <- X_var[, "EPU_State"]
  X_var_state <- X_var_state[complete.cases(X_var_state), ]
  X_var_state <- X_var_state[2:length(X_var_state)]
  X_var_state[X_var_state == 0] <-
    mean(X_var_state[X_var_state != 0])

  X_var_state <- makeReturns(X_var_state) * 100 * 12


  # M3: state-specific Composite EPU
  X_var_composite <- X_var[, "EPU_Composite"]
  X_var_composite <- X_var_composite[
    complete.cases(X_var_composite),
  ]
  X_var_composite <- X_var_composite[
    2:length(X_var_composite)
  ]
  X_var_composite[X_var_composite == 0] <-
    mean(X_var_composite[X_var_composite != 0])

  X_var_composite <-
    makeReturns(X_var_composite) * 100 * 12


  ############################################################
  # M4: average of the three state-specific EPU measures
  ############################################################

  X_var_all_mean <- apply(
    X_var,
    1,
    mean
  )

  X_var_all_mean <- as.xts(
    X_var_all_mean,
    time(X_var)
  )

  X_var_all_mean <- X_var_all_mean[
    complete.cases(X_var_all_mean),
  ]

  X_var_all_mean <- X_var_all_mean[
    2:length(X_var_all_mean)
  ]

  X_var_all_mean[X_var_all_mean == 0] <-
    mean(X_var_all_mean[X_var_all_mean != 0])

  X_var_all_mean <-
    makeReturns(X_var_all_mean) * 100 * 12

  ############################################################
  # M11: U.S. aggregate EPU
  ############################################################

  Full_X_var <- makeReturns(epu_full_i) * 100 * 12
  Full_X_var[1, ] <- 0

  ############################################################
  # M8-M10: Neighborhood EPU measures
  ############################################################

  # M8  : average National EPU of neighboring states
  # M9  : average Local EPU of neighboring states
  # M10 : average Composite EPU of neighboring states
  #
  # Alaska and Hawaii do not have contiguous U.S. neighbors.
  # For these two states, the corresponding state-specific EPU
  # measures are used instead.
  #
  # Maine has only one neighboring U.S. state, New Hampshire,
  # and is therefore handled separately.

  if (state_i %in% c("Alaska", "Hawaii")) {

    ##########################################################
    # Alaska and Hawaii
    ##########################################################

    Nat_neigh_mean_f   <- X_var_nat
    State_neigh_mean_f <- X_var_state
    Comp_neigh_mean_f  <- X_var_composite

  } else if (state_i == "Maine") {

    ##########################################################
    # Maine: New Hampshire is the only neighboring state
    ##########################################################

    X_var_neigh <- epu_by_state_i[["New Hampshire"]]

    # National EPU
    X_var_neigh_nat <- X_var_neigh[, "EPU_National"]
    X_var_neigh_nat <- X_var_neigh_nat[
      complete.cases(X_var_neigh_nat),
    ]
    X_var_neigh_nat <- X_var_neigh_nat[
      2:length(X_var_neigh_nat)
    ]
    X_var_neigh_nat[X_var_neigh_nat == 0] <-
      mean(X_var_neigh_nat[X_var_neigh_nat != 0])

    Nat_neigh_mean_f <-
      makeReturns(X_var_neigh_nat) * 100 * 12

    # Local EPU
    X_var_neigh_state <- X_var_neigh[, "EPU_State"]
    X_var_neigh_state <- X_var_neigh_state[
      complete.cases(X_var_neigh_state),
    ]
    X_var_neigh_state <- X_var_neigh_state[
      2:length(X_var_neigh_state)
    ]
    X_var_neigh_state[X_var_neigh_state == 0] <-
      mean(X_var_neigh_state[X_var_neigh_state != 0])

    State_neigh_mean_f <-
      makeReturns(X_var_neigh_state) * 100 * 12

    # Composite EPU
    X_var_neigh_composite <- X_var_neigh[, "EPU_Composite"]
    X_var_neigh_composite <- X_var_neigh_composite[
      complete.cases(X_var_neigh_composite),
    ]
    X_var_neigh_composite <- X_var_neigh_composite[
      2:length(X_var_neigh_composite)
    ]
    X_var_neigh_composite[X_var_neigh_composite == 0] <-
      mean(
        X_var_neigh_composite[
          X_var_neigh_composite != 0
        ]
      )

    Comp_neigh_mean_f <-
      makeReturns(X_var_neigh_composite) * 100 * 12

  } else {

    ##########################################################
    # All other contiguous states
    ##########################################################

    # Extract the adjacency vector for the current state.
    adj_mat_1 <- adj_mat[
      ordered_states_up[i],
    ]

    # Add Hawaii and Alaska, which are not included in usaww,
    # as zero-adjacency states.
    adj_mat_i <- c(
      adj_mat_1,
      HAWAII = 0,
      ALASKA = 0
    )

    # Reorder the adjacency vector according to the
    # population-based state ordering.
    adj_mat_i_f <- adj_mat_i[
      ordered_states_up
    ]

    # Identify neighboring states.
    vec_i <- which(
      adj_mat_i_f > 0
    )

    list_of_neighborhood <- epu_by_state_i[
      vec_i
    ]

    ##########################################################
    # Define the period needed for the MIDAS regressors
    ##########################################################

    # Five additional years of monthly data are retained before
    # the beginning of the daily return sample to accommodate the
    # MIDAS lag structure.
    year_b <- first(year(r_t)) - 5
    year_e <- 2024

    period_of_int <- paste0(
      year_b,
      "/",
      year_e
    )

    NN <- length(
      year_b:year_e
    ) * 12

    ##########################################################
    # Collect neighboring-state EPU series
    ##########################################################

    X_var_neigh_array <- array(
      NA,
      dim = c(
        NN,
        length(vec_i),
        3
      )
    )

    for (j in seq_along(vec_i)) {

      X_var_neigh <- list_of_neighborhood[[j]]

      # Preserve the original implementation by setting
      # the first observation to zero before extraction.
      X_var_neigh[1, ] <- 0

      X_var_neigh_array[, j, 1] <-
        X_var_neigh[, "EPU_National"][
          period_of_int
        ]

      X_var_neigh_array[, j, 2] <-
        X_var_neigh[, "EPU_State"][
          period_of_int
        ]

      X_var_neigh_array[, j, 3] <-
        X_var_neigh[, "EPU_Composite"][
          period_of_int
        ]
    }

    X_var_neigh_array[
      is.nan(X_var_neigh_array) |
      is.infinite(X_var_neigh_array)
    ] <- NA

    ##########################################################
    # Average across neighboring states
    ##########################################################

    # Use the dates of the state-level EPU series as reference.
    dates_neigh <- time(
      X_var[period_of_int]
    )

    if (length(vec_i) > 1) {

      Nat_neigh_mean <- apply(
        X_var_neigh_array[, , 1],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

      State_neigh_mean <- apply(
        X_var_neigh_array[, , 2],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

      Comp_neigh_mean <- apply(
        X_var_neigh_array[, , 3],
        1,
        function(x) mean(x, na.rm = TRUE)
      )

    } else {

      Nat_neigh_mean <-
        X_var_neigh_array[, , 1]

      State_neigh_mean <-
        X_var_neigh_array[, , 2]

      Comp_neigh_mean <-
        X_var_neigh_array[, , 3]
    }

    Nat_neigh_mean <- as.xts(
      Nat_neigh_mean,
      dates_neigh
    )

    State_neigh_mean <- as.xts(
      State_neigh_mean,
      dates_neigh
    )

    Comp_neigh_mean <- as.xts(
      Comp_neigh_mean,
      dates_neigh
    )

    ##########################################################
    # Convert to annualized percentage log-differences
    ##########################################################

    Nat_neigh_mean_f <-
      makeReturns(Nat_neigh_mean) * 100 * 12

    State_neigh_mean_f <-
      makeReturns(State_neigh_mean) * 100 * 12

    Comp_neigh_mean_f <-
      makeReturns(Comp_neigh_mean) * 100 * 12

    # Preserve the original implementation.
    Nat_neigh_mean_f[1, ]   <- 0
    State_neigh_mean_f[1, ] <- 0
    Comp_neigh_mean_f[1, ]  <- 0
  }

  ############################################################
  # Rolling-window estimation settings
  ############################################################

  # Number of daily return observations available for the
  # current state.
  TT <- length(r_t)

  # Number of rolling estimation steps.
  # Each model is estimated using a rolling window of Tin = 1000
  # daily observations and re-estimated every lstep = 25 days.
  nstep <- floor(
    (TT - Tin) / lstep
  )

  # Starting offsets of the rolling estimation windows.
  j <- c(
    0,
    (1:(nstep - 1)) * lstep
  )

  # First observation of each rolling estimation window.
  j_1 <- 1 + j

  # Last observation of each rolling estimation window.
  TTT <- Tin + j

  ############################################################
  # Initialize state-specific result containers
  ############################################################

  # The first column stores realized returns, while columns
  # 2-16 contain volatility estimates from M1-M15.
  in_sample_by_state_array <- array(
    NA,
    dim = c(
      Tin,
      16,
      nstep
    )
  )

  # Initial in-sample estimates.
  in_sample_by_state[[i]] <- matrix(
    NA,
    nrow = Tin,
    ncol = 16
  )

  in_sample_by_state[[i]][, 1] <- r_t[1:Tin]

  # Out-of-sample volatility forecasts.
  oos_by_state[[i]] <- matrix(
    NA,
    nrow = nstep * lstep,
    ncol = 16
  )

  # Multi-step-ahead volatility forecasts.
  oos_by_state_multi_step[[i]] <- matrix(
    NA,
    nrow = nstep * lstep,
    ncol = 16
  )

  # Store realized returns in the first column of both
  # out-of-sample result matrices.
  oos_by_state[[i]][, 1] <-
    r_t[(Tin + 1):(Tin + nstep * lstep)]

  oos_by_state_multi_step[[i]][, 1] <-
    r_t[(Tin + 1):(Tin + nstep * lstep)]


model_names <- c(
      "State-specific National EPU",
      "State-specific Local EPU",
      "State-specific Composite EPU",
      "Average state-specific EPU",
      "Cross-state average National EPU",
      "Cross-state average Local EPU",
      "Cross-state average Composite EPU",
      "Neighboring-state National EPU",
      "Neighboring-state Local EPU",
      "Neighboring-state Composite EPU",
      "U.S. aggregate EPU",
      "Realized Volatility",
      "Industrial Production",
      "Housing Starts"
    )

  ############################################################
  # Rolling in-sample and out-of-sample estimation loop
  ############################################################

  for (tt in 1:nstep) {

    ##########################################################
    # Define the current estimation and forecast windows
    ##########################################################

    # Last observation of the current in-sample estimation window.
    day_end_est <- TTT[tt]

    # First observation of the corresponding out-of-sample period.
    day_begin_oos <- day_end_est + 1

    # Last observation of the out-of-sample period.
    # For the final rolling step, use all remaining observations.
    if (tt == nstep) {

      day_end_oos <- TT

    } else {

      day_end_oos <- day_end_est + lstep
    }

    # First observation of the current rolling estimation window.
    day_begin_est <- j_1[tt]

    ##########################################################
    # Construct return samples for the current rolling step
    ##########################################################

    # Return series used by ugmfit().
    #
    # The object includes the Tin in-sample observations plus
    # the following lstep observations required for the
    # out-of-sample forecasts.
    r_t_est_cycle <- r_t[
      day_begin_est:(day_end_est + lstep)
    ]

    # Current rolling in-sample window.
    r_t_in_s <- r_t[
      day_begin_est:day_end_est
    ]

    # Last in-sample return, used for multi-step GJR forecasts.
    r_t_gjr_multi <- last(r_t_in_s)

    # Returns corresponding to the next lstep forecast observations.
    r_t_oos <- r_t[
      day_begin_oos:(day_begin_oos + lstep - 1)
    ]

    # Return series used to align the GARCH-MIDAS variables.
    # Unlike r_t_est_cycle, this series always starts from the
    # beginning of the state-specific sample.
    r_t_est_cycle_gm <- r_t[
      1:(day_end_est + lstep)
    ]

    ##########################################################
    # Store current in-sample returns
    ##########################################################

    in_sample_by_state_array[, 1, tt] <-
      r_t_in_s


    ##########################################################
    # MIDAS variables
    ##########################################################

    # Construct the MIDAS matrices for the 14 GARCH-MIDAS
    # specifications considered in the forecasting exercise.

    # M1-M4: state-specific EPU measures
    m1_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      X_var_nat,
      K = K,
      "monthly"
    )

    m2_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      X_var_state,
      K = K,
      "monthly"
    )

    m3_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      X_var_composite,
      K = K,
      "monthly"
    )

    m4_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      X_var_all_mean,
      K = K,
      "monthly"
    )

    # M5-M7: cross-state average EPU measures
    m5_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      Nat_mean,
      K = K,
      "monthly"
    )

    m6_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      State_mean,
      K = K,
      "monthly"
    )

    m7_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      Comp_mean,
      K = K,
      "monthly"
    )

    # M8-M10: neighboring-state average EPU measures
    m8_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      Nat_neigh_mean_f,
      K = K,
      "monthly"
    )

    m9_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      State_neigh_mean_f,
      K = K,
      "monthly"
    )

    m10_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      Comp_neigh_mean_f,
      K = K,
      "monthly"
    )

    # M11: U.S. aggregate EPU
    m11_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      Full_X_var,
      K = K,
      "monthly"
    )

    # M12-M14: alternative macro-financial predictors
    m12_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      RV_monthly,
      K = K,
      "monthly"
    )

    m13_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      ip_f_i,
      K = K,
      "monthly"
    )

    m14_mv_cycle <- mv_into_mat(
      r_t_est_cycle,
      hs_f_i,
      K = K,
      "monthly"
    )

    ##########################################################
    # Collect MIDAS matrices for M1-M14
    ##########################################################

    midas_variables <- list(
      m1_mv_cycle,
      m2_mv_cycle,
      m3_mv_cycle,
      m4_mv_cycle,
      m5_mv_cycle,
      m6_mv_cycle,
      m7_mv_cycle,
      m8_mv_cycle,
      m9_mv_cycle,
      m10_mv_cycle,
      m11_mv_cycle,
      m12_mv_cycle,
      m13_mv_cycle,
      m14_mv_cycle
    )

    
# m15: gjr

    ##########################################################
    # Estimate M1-M14: GARCH-MIDAS models
    ##########################################################

    for (m in 1:14) {

      # Columns 2-15 of the result matrices correspond to M1-M14.
      COL <- m + 1

      message(
        "State: ", state_i,
        " | rolling step: ", tt, "/", nstep,
        " | model M", m, ": ", model_names[m]
      )

      ########################################################
      # Estimate GARCH-MIDAS model
      ########################################################

      count <- 1

      repeat {

        gm_est <- NULL

        gm_est <- tryCatch(
          ugmfit(
            model = "GM",
            skew = "YES",
            distribution = "norm",
            daily_ret = r_t_est_cycle,
            mv_m = midas_variables[[m]],
            K = K,
            R = 1000,
            out_of_sample = lstep
          ),
          error = function(e) return(NA)
        )

        count <- count + 1

        # Stop when estimation succeeds or after the maximum
        # number of attempts used in the original implementation.
        if (all(!is.na(gm_est)) | count == 10) {
          break
        }
      }

      ########################################################
      # Store results from the initial estimation window
      ########################################################

      if (tt == 1) {

        est_coef_mat <- gm_est$rob_coef_mat

        # Statistical significance of the MIDAS slope coefficient.
        sig <- ifelse(
          est_coef_mat[5, 4] <= 0.01,
          "^{***}",
          ifelse(
            est_coef_mat[5, 4] <= 0.05,
            "^{**}",
            ifelse(
              est_coef_mat[5, 4] <= 0.10,
              "^{*}",
              "^{}"
            )
          )
        )

        # Store the MIDAS slope coefficient and significance level.
        tab_mat_est[i, COL] <- paste0(
          round(est_coef_mat[5, 1], 3),
          sig
        )

        # Store estimated in-sample volatility from the
        # initial rolling window.
        in_sample_by_state[[i]][, COL] <-
          coredata(gm_est$est_vol_in_s)
      }

      ########################################################
      # Define positions of the current OOS forecasts
      ########################################################

      if (tt < nstep) {

        oos_index <- j_1[tt]:(j_1[tt + 1] - 1)

      } else {

        oos_index <- j_1[tt]:(nstep * lstep)
      }

      ########################################################
      # Store one-step-ahead forecasts
      ########################################################

      oos_by_state[[i]][oos_index, COL] <-
        coredata(gm_est$est_vol_oos)

      ########################################################
      # Store multi-step-ahead forecasts
      ########################################################

      oos_by_state_multi_step[[i]][oos_index, COL] <-
        multi_step_ahead_pred(
          gm_est,
          lstep
        )

      ########################################################
      # Store rolling-window in-sample volatility estimates
      ########################################################

      in_sample_by_state_array[, COL, tt] <-
        coredata(gm_est$est_vol_in_s)
    }

    ##########################################################
    # M15: GJR-GARCH 
    ##########################################################

    # Column 16 corresponds to the GJR-GARCH.
    COL <- 16

    message(
      "State: ", state_i,
      " | rolling step: ", tt, "/", nstep,
      " | model M15: GJR-GARCH"
    )

    ##########################################################
    # Estimate GJR-GARCH
    ##########################################################

    fit_gjr_n <- ugarchfit(
      spec = spec_gjr_n,
      data = r_t_est_cycle,
      out.sample = lstep,
      solver = "hybrid"
    )

    ##########################################################
    # In-sample volatility
    ##########################################################

    sigma_gjr_n_in_s <-
      fit_gjr_n@fit$sigma

    ##########################################################
    # One-step-ahead rolling forecasts
    ##########################################################

    sigma_gjr_n_oos <- as.numeric(
      sigma(
        ugarchforecast(
          fit_gjr_n,
          n.ahead = 1,
          n.roll = lstep - 1,
          out.sample = lstep - 1
        )
      )
    )

    ##########################################################
    # Store initial in-sample estimates
    ##########################################################

    if (tt == 1) {

      in_sample_by_state[[i]][, COL] <-
        sigma_gjr_n_in_s
    }

    ##########################################################
    # Store out-of-sample forecasts
    ##########################################################

    # oos_index is the same index used above for M1-M14.
    oos_by_state[[i]][oos_index, COL] <-
      sigma_gjr_n_oos

    # Multi-step-ahead GJR-GARCH forecasts.
    oos_by_state_multi_step[[i]][oos_index, COL] <-
      gjr_multi_step_ahead(
        fit_gjr_n,
        r_t_gjr_multi,
        lstep
      )

    ##########################################################
    # Store rolling-window in-sample volatility estimates
    ##########################################################

    in_sample_by_state_array[, COL, tt] <-
      sigma_gjr_n_in_s

    ##########################################################
    # Save intermediate results
    ##########################################################

    # Save the current rolling-step results so that intermediate
    # estimates are available if the full estimation is interrupted.
    save(
      tt,
      oos_by_state,
      oos_by_state_multi_step,
      file = temp_filename
    )

    ##########################################################
    # End of rolling in-sample / out-of-sample cycle
    ##########################################################

    cat(
      "State:", state_i,
      "| step", tt,
      "out of", nstep,
      "\n"
    )
  }

  ############################################################
  # Store state-specific rolling in-sample results
  ############################################################

  in_sample_by_state_list_array[[i]] <-
    in_sample_by_state_array

  message(
    "Completed state ", i,
    " of 50: ", state_i
  )

  ############################################################
  # Save final results for the current state subset
  ############################################################

  save(
    i,
    tt,
    in_sample_by_state,
    oos_by_state,
    oos_by_state_multi_step,
    tab_mat_est,
    in_sample_by_state_list_array,
    file = filename
  )

}  # end state-level estimation loop

####################################### 
####################################### END cycle
#######################################


############################################################
# Elastic Net forecast combination
############################################################

# The forecast combination is constructed using the volatility
# forecasts from the EPU-based GARCH-MIDAS specifications
# (M1-M11) together with the GJR-GARCH benchmark.
#
# The following GARCH-MIDAS specifications are therefore excluded:
#   M12: Realized Volatility
#   M13: Industrial Production
#   M14: Housing Starts
#
# In the result matrices:
#   column 1     = observed return
#   columns 2-12 = EPU-based GARCH-MIDAS models M1-M11
#   columns 13-15 = RV, IP, and HS GARCH-MIDAS models
#   column 16    = GJR-GARCH
#
# Hence, the Elastic Net combination uses columns 2:12 and 16.
#
# Absolute returns are used as a proxy for latent volatility.
# The individual-model volatility estimates and forecasts are
# used in their original annualized scale, with no deannualization.

############################################################
# Combination settings
############################################################

# Elastic Net mixing parameter.
alpha_f <- 0.5


############################################################
# Lambda grids
############################################################

lambda_grid <- expand.grid(
  lambda = 20^seq(-4, 1, length = 100),
  alpha = alpha_f
)

# Alternative grid used when the lambda selected from the
# primary grid is located at its upper boundary.
lambda_grid_2 <- expand.grid(
  lambda = 10^seq(-4, 1, length = 100),
  alpha = alpha_f
)


############################################################
# Time-series validation scheme
############################################################

# Lambda is selected using rolling time-series validation.
# For each 1000-observation in-sample window:
#   - the initial training window contains 800 observations;
#   - validation is performed one observation ahead at a time;
#   - a fixed-width training window is used.

train_control <- trainControl(
  method = "timeslice",
  initialWindow = 800,
  horizon = 1,
  fixedWindow = TRUE,
  skip = 0,
  verboseIter = FALSE
)


############################################################
# Initialize result containers
############################################################

# Copies of the original OOS forecast objects. The Elastic Net
# combination forecasts are appended as additional columns.

oos_by_state_2 <- oos_by_state

oos_by_state_multi_step_2 <-
  oos_by_state_multi_step

# Selected lambda values and estimated Elastic Net coefficients.

lambda_list <- list()
coef_est <- list()

# EPU-based predictor set

# M1-M11 (columns 2:12) + GJR-GARCH (column 16).
epu_model_cols <- c(2:12, 16)


############################################################
# State-level forecast-combination cycle
############################################################

for (i in i_b:i_e) {

  set.seed(i)

  state_i <- ordered_states[i]

  message(
    "Estimating Elastic Net combination for state ",
    i, " of 50: ", state_i
  )

  ############################################################
  # Extract state-specific in-sample and OOS forecasts
  ############################################################

  # The in-sample array contains:
  #   column 1     = observed returns
  #   columns 2-15 = GARCH-MIDAS models M1-M14
  #   column 16    = GJR-GARCH
  #
  # The Elastic Net combination uses M1-M11 and GJR-GARCH.
  # Hence, columns 13-15 (RV, IP, and HS) are excluded.

  tab_in_s <- in_sample_by_state_list_array[[i]][
  	, c(1, epu_model_cols),
  ]

  # For the OOS objects, column 1 contains observed returns and
  # is not required for constructing the combination forecasts.

  tab_oos <- oos_by_state[[i]][
    , epu_model_cols
  ]

  tab_oos_ms <- oos_by_state_multi_step[[i]][
    , epu_model_cols
  ]


  ############################################################
  # Rolling-window settings
  ############################################################

  # Number of rolling estimation windows available for
  # the current state.
  nstep <- dim(
     in_sample_by_state_list_array[[i]]
  )[3]

  # Starting positions of the 25-observation OOS blocks.
  j <- c(
    0,
    (1:(nstep - 1)) * lstep
  )

  j_1 <- 1 + j


  ############################################################
  # Initialize state-specific result containers
  ############################################################

  # Selected lambda for each rolling window.
  lambda_l <- list()

  # Elastic Net coefficient estimates for each rolling window.
  coef_e <- list()

  # One-step-ahead and multi-step-ahead combination forecasts.
  y_hat_one_step <- list()
  y_hat_multi_step <- list()


  ############################################################
  # Rolling Elastic Net combination cycle
  ############################################################

  for (tt in 1:nstep) {

    ##########################################################
    # Define current OOS block
    ##########################################################

    if (tt < nstep) {

      oos_index <- j_1[tt]:(j_1[tt + 1] - 1)

    } else {

      oos_index <- j_1[tt]:nrow(tab_oos)
    }


    ##########################################################
    # In-sample and OOS predictor matrices
    ##########################################################

    # Current 1000-observation in-sample window.
    #
    # Column 1 contains observed returns, while the remaining
    # columns contain the volatility estimates from M1-M11
    # and GJR-GARCH.

    X_in_s <- tab_in_s[
      , , tt
    ]

    # Individual-model OOS volatility forecasts used to
    # construct the one-step-ahead combination.

    X_oos <- tab_oos[
      oos_index,
    ]

    # Individual-model multi-step-ahead volatility forecasts.

    X_oos_ms <- tab_oos_ms[
      oos_index,
    ]


    ##########################################################
    # Construct Elastic Net target and predictors
    ##########################################################

    # Predictors: in-sample volatility estimates from
    # M1-M11 and GJR-GARCH.

    X <- X_in_s[
      , 2:ncol(X_in_s)
    ]

    # Target: absolute returns, used as a proxy for
    # the latent daily volatility.

    y <- abs(
      X_in_s[, 1]
    )

    # Assign simple predictor names for caret/glmnet.
    colnames(X) <- seq_len(
      ncol(X)
    )

    y <- as.numeric(
      coredata(y)
    )


    ##########################################################
    # Select lambda using time-series validation
    ##########################################################

    fit <- train(
      x = X,
      y = y,
      method = "glmnet",
      tuneGrid = lambda_grid,
      trControl = train_control,
      lower.limits = 0
    )

    best_lambda <- fit$bestTune$lambda

    # If the selected lambda is the upper boundary of the
    # primary grid, repeat the search using the alternative grid.

    if (best_lambda == 20) {

      fit <- train(
        x = X,
        y = y,
        method = "glmnet",
        tuneGrid = lambda_grid_2,
        trControl = train_control,
        lower.limits = 0
      )

      best_lambda <- fit$bestTune$lambda
    }

    lambda_l[[tt]] <- best_lambda


    ##########################################################
    # Estimate Elastic Net using the selected lambda
    ##########################################################

    # Non-negative constraints are imposed on the coefficients
    # associated with the individual volatility forecasts.

    enet_model_best <- glmnet(
      X,
      y,
      alpha = alpha_f,
      lambda = best_lambda,
      intercept = TRUE,
      lower.limits = 0
    )

    coef_e[[tt]] <- coef(
      enet_model_best
    )


    ##########################################################
    # Construct Elastic Net combination forecasts
    ##########################################################

    y_hat_one_step[[tt]] <- predict(
      enet_model_best,
      X_oos
    )

    y_hat_multi_step[[tt]] <- predict(
      enet_model_best,
      X_oos_ms
    )


    ##########################################################
    # Progress message
    ##########################################################

    message(
      "State: ", state_i,
      " | combination step: ",
      tt, "/", nstep
    )
  }


  ############################################################
  # Store state-specific Elastic Net results
  ############################################################

  lambda_list[[i]] <- lambda_l

  coef_est[[i]] <- coef_e


  ############################################################
  # Append Elastic Net forecasts to OOS result matrices
  ############################################################

  oos_by_state_2[[i]] <- cbind(
    oos_by_state[[i]],
    EN_Comb = unlist(
      y_hat_one_step
    )
  )

  oos_by_state_multi_step_2[[i]] <- cbind(
    oos_by_state_multi_step[[i]],
    EN_Comb = unlist(
      y_hat_multi_step
    )
  )


  ############################################################
  # Save intermediate results
  ############################################################

  filename_comb <- paste0(
    "comb_lstep_25_", i_b, "_", i_e, ".RData"
  )

  save(
    i,
    lambda_list,
    coef_est,
    oos_by_state_2,
    oos_by_state_multi_step_2,
    file = filename_comb
  )

  message(
    "Completed Elastic Net combination for ",
    state_i
  )

}  # end state-level combination cycle



############################################################
# Table 3: Percentage of non-zero Elastic Net coefficients
############################################################

# Table 3 reports, for each state and each predictor included in
# the Elastic Net forecast combination, the percentage of rolling
# estimation windows in which the corresponding coefficient is
# different from zero.
#
# The Elastic Net combination includes:
#   - M1-M11: EPU-based GARCH-MIDAS specifications;
#   - M15: GJR-GARCH benchmark.
#
# The intercept is excluded when computing selection frequencies.

lab_tab <- lab <- c(
  "GARCH-MIDAS (Nat.)",
  "GARCH-MIDAS (Loc.)",
  "GARCH-MIDAS (Comp.)",
  "GARCH-MIDAS (Avg.)",
  "GARCH-MIDAS (Nat. Avg.)",
  "GARCH-MIDAS (Loc. Avg.)",
  "GARCH-MIDAS (Comp. Avg.)",
  "GARCH-MIDAS (Neigh. Nat.)",
  "GARCH-MIDAS (Neigh. Loc.)",
  "GARCH-MIDAS (Neigh. Comp.)",
  "GARCH-MIDAS (Glob.)",
  "GJR"
)

############################################################
# Initialize selection-frequency matrix
############################################################

# One row for each U.S. state and one column for each predictor
# included in the Elastic Net combination.
perc_selected <- matrix(
  NA,
  nrow = 50,
  ncol = 12
)

colnames(perc_selected) <- lab_tab

# Number of rolling estimation windows available for each state.
obs <- rep(
  NA,
  50
)


############################################################
# Compute state-specific selection frequencies
############################################################

for (i in i_b:i_e) {

  # Elastic Net coefficient estimates across rolling windows
  # for state i.
  X_state <- coef_est[[i]]

  # Number of rolling estimation windows for the current state.
  nstep <- obs[i] <- length(X_state)

  # For each rolling window, identify whether each coefficient
  # is different from zero.
  #
  # The first coefficient is the intercept and is therefore removed.
  sel_mat <- sapply(
    1:nstep,
    function(tt) {

      beta <- X_state[[tt]][-1]

      beta != 0
    }
  )

  # Percentage of rolling windows in which each predictor
  # receives a non-zero Elastic Net coefficient.
  perc_selected[i, ] <- rowMeans(
    sel_mat
  ) * 100
}


############################################################
# Cross-state averages
############################################################

# Weighted cross-state average reported in the table.
#
# States are weighted by the number of rolling estimation windows
# available for the corresponding state. This accounts for the
# different sample lengths across states.
w_mean <- 100 *
  colSums(
    (perc_selected / 100) * obs,
    na.rm = TRUE
  ) /
  sum(
    obs,
    na.rm = TRUE
  )


############################################################
# Construct Table 3
############################################################

# Add the weighted cross-state average as the final row.
tab_matrix <- rbind(
  perc_selected,
  w_mean
)

# State names follow the population-based ordering used throughout
# the empirical analysis.
row_names <- c(
  ordered_states,
  "Average"
)

# Construct the final table.
table_3 <- data.frame(
  State = row_names,
  Obs = c(obs, NA),
  formatC(
    tab_matrix,
    format = "f",
    digits = 3
  ),
  check.names = FALSE
)

table_3


############################################################
# LaTeX output
############################################################

xtable(
  table_3,
  type = "latex"
) |>
  print(
    include.rownames = FALSE
  )


############################################################
# Assign model names to OOS forecast matrices
############################################################

model_colnames <- c(
  "Return",
  "M1_Nat",
  "M2_Loc",
  "M3_Comp",
  "M4_Avg",
  "M5_Mean_Nat",
  "M6_Mean_Loc",
  "M7_Mean_Comp",
  "M8_Neigh_Nat",
  "M9_Neigh_Loc",
  "M10_Neigh_Comp",
  "M11_US_EPU",
  "M12_RV",
  "M13_IP",
  "M14_HS",
  "M15_GJR",
  "EN_Comb"
)

for (i in i_b:i_e) {

  colnames(oos_by_state_2[[i]]) <-
    model_colnames

  colnames(oos_by_state_multi_step_2[[i]]) <-
    model_colnames
}


############################################################
# Mean and Median forecast combinations
############################################################

# The Mean and Median combinations are computed using all
# 15 individual volatility models (M1-M14 and GJR-GARCH).
#
# The Elastic Net combination is not included among the
# individual forecasts entering these benchmark combinations.

individual_models <- model_colnames[2:16]


############################################################
# State-level combination cycle
############################################################

for (i in i_b:i_e) {

  state_i <- ordered_states[i]

  ############################################################
  # Individual-model OOS forecasts
  ############################################################

  # One-step-ahead forecasts.
  X_ss <- oos_by_state_2[[i]][
    , individual_models
  ]

  # Multi-step-ahead forecasts.
  X_ms <- oos_by_state_multi_step_2[[i]][
    , individual_models
  ]


  ############################################################
  # Mean combination
  ############################################################

  mean_comb_ss <- apply(
    X_ss,
    1,
    mean,
    na.rm = TRUE
  )

  mean_comb_ms <- apply(
    X_ms,
    1,
    mean,
    na.rm = TRUE
  )


  ############################################################
  # Median combination
  ############################################################

  median_comb_ss <- apply(
    X_ss,
    1,
    median,
    na.rm = TRUE
  )

  median_comb_ms <- apply(
    X_ms,
    1,
    median,
    na.rm = TRUE
  )


  ############################################################
  # Append combination forecasts
  ############################################################

  oos_by_state_2[[i]] <- cbind(
    oos_by_state_2[[i]],
    Mean_Comb = mean_comb_ss,
    Median_Comb = median_comb_ss
  )

  oos_by_state_multi_step_2[[i]] <- cbind(
    oos_by_state_multi_step_2[[i]],
    Mean_Comb = mean_comb_ms,
    Median_Comb = median_comb_ms
  )


  ############################################################
  # Progress message
  ############################################################

  message(
    "Mean and Median combinations completed for ",
    state_i
  )
}

############################################################
# MSE- and QLIKE-based forecast combinations
############################################################

# The MSE and QLIKE combinations are constructed using all
# 15 individual volatility models:
#
#   M1-M14 + GJR-GARCH.
#
# The Elastic Net, Mean, and Median combinations are not
# included among the forecasts used to compute these weights.
#
# At each rolling estimation step, model-specific weights are
# computed from the corresponding 1000-observation in-sample
# window.
#
# Absolute returns are used as the volatility proxy.

############################################################
# Individual model set
############################################################

individual_models <- model_colnames[2:16]


############################################################
# State-level combination cycle
############################################################

for (i in i_b:i_e) {

  state_i <- ordered_states[i]

  message(
    "Computing MSE and QLIKE combinations for state ",
    i, " of 50: ", state_i
  )

  ############################################################
  # State-specific in-sample and OOS forecasts
  ############################################################

  # Rolling in-sample arrays:
  # column 1 contains returns;
  # columns 2-16 contain M1-M15 volatility estimates.
  tab_in_s <- in_sample_by_state_list_array[[i]]

  # Individual-model one-step-ahead OOS forecasts.
  tab_oos <- oos_by_state_2[[i]][
    , individual_models
  ]

  # Individual-model multi-step-ahead OOS forecasts.
  tab_oos_ms <- oos_by_state_multi_step_2[[i]][
    , individual_models
  ]

  ############################################################
  # Rolling-window settings
  ############################################################

  nstep <- dim(
    tab_in_s
  )[3]

  j <- c(
    0,
    (1:(nstep - 1)) * lstep
  )

  j_1 <- 1 + j

  ############################################################
  # Initialize state-specific combination forecasts
  ############################################################

  Comb_mse_ss_l <- list()
  Comb_mse_ms_l <- list()

  Comb_qlike_ss_l <- list()
  Comb_qlike_ms_l <- list()


  ############################################################
  # Rolling weighting cycle
  ############################################################

  for (tt in 1:nstep) {

    ##########################################################
    # Current in-sample window
    ##########################################################

    X_in_s <- tab_in_s[
      , , tt
    ]

    # Absolute returns are used as the volatility proxy.
    abs_ret <- abs(
      X_in_s[, 1]
    )

    # In-sample volatility estimates from M1-M15.
    all_models <- X_in_s[
      , 2:16
    ]


    ##########################################################
    # Compute model-specific MSE and QLIKE
    ##########################################################

    mse_tt <- colMeans(
      MSE_f(
        all_models,
        abs_ret
      )
    )

    qlike_tt <- colMeans(
      QLIKE_f(
        abs_ret,
        all_models
      )
    )


    ##########################################################
    # Construct inverse-loss weights
    ##########################################################

    inv_mse_tt <- mse_tt^(-1)
    inv_qlike_tt <- qlike_tt^(-1)

    # Normalize weights so that they sum to one.
    w_mse <- inv_mse_tt / sum(
      inv_mse_tt
    )

    w_qlike <- inv_qlike_tt / sum(
      inv_qlike_tt
    )


    ##########################################################
    # Define current OOS block
    ##########################################################

    if (tt < nstep) {

      oos_index <- j_1[tt]:(j_1[tt + 1] - 1)

    } else {

      oos_index <- j_1[tt]:nrow(tab_oos)
    }


    ##########################################################
    # MSE-weighted combinations
    ##########################################################

    Comb_mse_ss_l[[tt]] <- rowSums(
      tab_oos[
        oos_index,
      ] * w_mse
    )

    Comb_mse_ms_l[[tt]] <- rowSums(
      tab_oos_ms[
        oos_index,
      ] * w_mse
    )


    ##########################################################
    # QLIKE-weighted combinations
    ##########################################################

    Comb_qlike_ss_l[[tt]] <- rowSums(
      tab_oos[
        oos_index,
      ] * w_qlike
    )

    Comb_qlike_ms_l[[tt]] <- rowSums(
      tab_oos_ms[
        oos_index,
      ] * w_qlike
    )


    ##########################################################
    # Progress message
    ##########################################################

    message(
      "State: ", state_i,
      " | weighting step: ",
      tt, "/", nstep
    )
  }


  ############################################################
  # Append MSE and QLIKE combinations
  ############################################################

  oos_by_state_2[[i]] <- cbind(
    oos_by_state_2[[i]],
    MSE_Comb = unlist(
      Comb_mse_ss_l
    ),
    QLIKE_Comb = unlist(
      Comb_qlike_ss_l
    )
  )

  oos_by_state_multi_step_2[[i]] <- cbind(
    oos_by_state_multi_step_2[[i]],
    MSE_Comb = unlist(
      Comb_mse_ms_l
    ),
    QLIKE_Comb = unlist(
      Comb_qlike_ms_l
    )
  )

  message(
    "MSE and QLIKE combinations completed for ",
    state_i
  )

}  # end state-level combination cycle


############################################################
# Tables 4-7: Out-of-sample forecast evaluation
############################################################

# Four tables are constructed:
#
# Table 4: QLIKE loss, one-step-ahead forecasts
# Table 5: MSE loss,   one-step-ahead forecasts
# Table 6: QLIKE loss, multi-step-ahead forecasts
# Table 7: MSE loss,   multi-step-ahead forecasts
#
# All returns and volatility forecasts are converted back from
# annualized percentage units to their original daily scale before
# computing the loss functions.
#
# QLIKE values are reported in their original scale.
# MSE values are multiplied by 10,000 to improve readability.
#
# For each state, the Model Confidence Set (MCS) is computed using
# the squared-statistic version returned by MCS_f().

alpha <- 0.25 	# significance level of the MCS
B <- 5000		# number of bootstrap replicated used in MCS

# Convert annualized percentage returns/volatilities back to
# daily non-percentage units.
coeff_deann <- (100 * sqrt(252))^(-1)


############################################################
# Model labels
############################################################

model_labels <- c(
  "Nat.",
  "State",
  "Comp.",
  "Mean",
  "Mean-Nat.",
  "Mean-State",
  "Mean-Comp.",
  "Mean-Nat.-Neigh.",
  "Mean-State-Neigh.",
  "Mean-Comp.-Neigh.",
  "Full",
  "RV",
  "IP",
  "HS",
  "GJR",
  "EN-Comb",
  "Mean-Comb",
  "Median-Comb",
  "MSE-Comb",
  "QLIKE-Comb"
)


############################################################
# Table specifications
############################################################

table_specs <- list(

  table_4 = list(
    forecasts = oos_by_state_2,
    loss = QLIKE_f,
    coeff = 1,
    description = "QLIKE, one-step-ahead"
  ),

  table_5 = list(
    forecasts = oos_by_state_2,
    loss = MSE_f,
    coeff = 10000,
    description = "MSE, one-step-ahead"
  ),

  table_6 = list(
    forecasts = oos_by_state_multi_step_2,
    loss = QLIKE_f,
    coeff = 1,
    description = "QLIKE, multi-step-ahead"
  ),

  table_7 = list(
    forecasts = oos_by_state_multi_step_2,
    loss = MSE_f,
    coeff = 10000,
    description = "MSE, multi-step-ahead"
  )
)


############################################################
# Construct Tables 4-7
############################################################

results_tables <- list()

for (table_name in names(table_specs)) {

  spec <- table_specs[[table_name]]

  forecast_object <- spec$forecasts
  LOSS <- spec$loss
  COEFF <- spec$coeff

  message(
    "Constructing ", table_name,
    " (", spec$description, ")"
  )


  ##########################################################
  # Initialize state-level table
  ##########################################################

  NC <- ncol(
    forecast_object[[1]]
  )

  col_mcs <- as.data.frame(
    matrix(
      NA,
      nrow = 50,
      ncol = NC
    )
  )

  col_mcs[, 1] <- ordered_states

  colnames(col_mcs) <- c(
    "State",
    model_labels
  )


  ##########################################################
  # State-level MCS cycle
  ##########################################################

  for (j in i_b:i_e) {

    state_i <- ordered_states[j]

    ########################################################
    # Prepare returns and forecasts
    ########################################################

    # Take absolute values and convert all series back to
    # daily non-percentage units.
    tab_f <- abs(
      coredata(
        forecast_object[[j]]
      )
    ) * coeff_deann

    # Column 1 contains absolute returns.
    # Columns 2-21 contain the 20 volatility forecasts.
    target <- tab_f[, 1]

    forecasts <- tab_f[
      , 2:ncol(tab_f),
      drop = FALSE
    ]


    ########################################################
    # Handle missing forecast values
    ########################################################

    # If one or more model forecasts are missing for a given
    # observation, replace the missing values with the mean
    # of the available forecasts for that observation.
    forecasts <- t(
      apply(
        forecasts,
        1,
        function(x) {

          if (any(is.na(x))) {

            row_mean <- mean(
              x,
              na.rm = TRUE
            )

            x[is.na(x)] <- row_mean
          }

          x
        }
      )
    )

    ########################################################
    # Compute model-specific losses
    ########################################################

    db_mcs <- sapply(
      1:ncol(forecasts),
      function(k) {

        LOSS(
          target,
          forecasts[, k]
        )
      }
    )

    db_mcs <- as.matrix(
      db_mcs
    )

    colnames(db_mcs) <- model_labels


    ########################################################
    # Model Confidence Set
    ########################################################


	table_id <- match(
  	  table_name,
      names(table_specs)
    )

	  set.seed(
    1000 * table_id + j
    )

    mcs_res <- MCS_f(
      db_mcs,
      B,
      alpha
    )

    # Use the MCS based on the squared statistic.
    mcs_to_use <- mcs_res$includedSQ


    ########################################################
    # Average loss and MCS membership
    ########################################################

    col_means <- colMeans(
      db_mcs
    )

    col_mcs_2 <- cbind(
      round(
        col_means * COEFF,
        3
      ),
      ifelse(
        1:ncol(db_mcs) %in% mcs_to_use,
        1,
        0
      )
    )

    # Prefix values belonging to the MCS with "cell".
    # This marker is subsequently used by DT for highlighting.
    row_mcs_f <- paste0(
      ifelse(
        col_mcs_2[, 2] == 1,
        "cell",
        ""
      ),
      col_mcs_2[, 1]
    )

    col_mcs[
      j,
      2:ncol(col_mcs)
    ] <- row_mcs_f


    message(
      table_name,
      " | state ",
      j, "/50: ",
      state_i
    )
  }


  ##########################################################
  # Store completed table
  ##########################################################

  results_tables[[table_name]] <- col_mcs
}


############################################################
# Create the four final table objects
############################################################

table_4 <- results_tables$table_4
table_5 <- results_tables$table_5
table_6 <- results_tables$table_6
table_7 <- results_tables$table_7

############################################################
# HTML visualization
############################################################

display_mcs_table(table_4)
display_mcs_table(table_5)
display_mcs_table(table_6)
display_mcs_table(table_7)

