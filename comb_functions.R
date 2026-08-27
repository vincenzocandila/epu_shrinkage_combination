############################################################
# Summary statistics
############################################################

summ_stat_f <- function(x) {

  # Computes basic descriptive statistics for a numeric vector.
  #
  # Input:
  #   x : numeric vector.
  #
  # Output:
  #   A one-row matrix containing:
  #     col_1 : minimum;
  #     col_2 : maximum;
  #     col_3 : mean;
  #     col_4 : standard deviation.

  col_1 <- min(x)
  col_3 <- mean(x)
  col_4 <- sd(x)
  col_2 <- max(x)

  return(
    cbind(
      col_1,
      col_2,
      col_3,
      col_4
    )
  )
}


############################################################
# DT visualization
############################################################

display_mcs_table <- function(x) {

  # Displays an MCS result table using DT.
  #
  # Input:
  #   x : data frame containing the results of the forecast
  #       evaluation. Entries prefixed by "cell" identify models
  #       belonging to the Model Confidence Set.
  #
  # Output:
  #   A DT datatable in which entries belonging to the MCS
  #   are highlighted in light green.

  x_styled <- as.data.frame(
    lapply(
      x,
      function(z) {

        ifelse(
          grepl("^cell", z),
          paste0(
            '<div style="background-color: lightgreen;">',
            gsub("^cell", "", z),
            "</div>"
          ),
          gsub("^cell", "", z)
        )
      }
    ),
    stringsAsFactors = FALSE
  )

  datatable(
    x_styled,
    escape = FALSE,
    options = list(
      scrollX = FALSE,
      paging = FALSE
    )
  )
}


############################################################
# GJR-GARCH multi-step-ahead forecasts
############################################################

gjr_multi_step_ahead <- function(x, RET, k) {

  # Computes multi-step-ahead volatility forecasts from an
  # estimated GJR-GARCH model.
  #
  # Inputs:
  #   x   : fitted GJR-GARCH model object returned by rugarch;
  #   RET : most recent return used to initialize the forecast;
  #   k   : forecast horizon.
  #
  # Output:
  #   Numeric vector of length k containing volatility forecasts
  #   for horizons 1,...,k.

  h_t <- last(x@fit$sigma)^2

  alpha_0 <- coef(x)[1]
  alpha_1 <- coef(x)[2]
  beta_1  <- coef(x)[3]
  gamma_1 <- coef(x)[4]

  res <- numeric(k)

  for (j in 1:k) {

    term1 <- sum(
      alpha_0 *
        ((alpha_1 + gamma_1 / 2 + beta_1) ^ (0:(j - 1)))
    )

    term2 <-
      (alpha_1 + gamma_1 / 2 + beta_1)^(j - 1) *
      (
        alpha_1 * RET^2 +
        beta_1 * h_t +
        gamma_1 * RET^2 * (RET < 0)
      )

    res[j] <- term1 + term2
  }

  res_vol <- sqrt(res)

  return(res_vol)
}


############################################################
# QLIKE loss function
############################################################

QLIKE_f <- function(x, y) {

  # Computes the QLIKE loss.
  #
  # Inputs:
  #   x : volatility proxy;
  #   y : volatility estimate or forecast.
  #
  # Output:
  #   QLIKE loss values with the same conformable dimensions
  #   as x and y.

  log(y) + (x / y)
}


############################################################
# MSE loss function
############################################################

MSE_f <- function(x, y) {

  # Computes the Mean Squared Error loss.
  #
  # Inputs:
  #   x : volatility proxy;
  #   y : volatility estimate or forecast.
  #
  # Output:
  #   Squared-error loss values with the same conformable
  #   dimensions as x and y.

  (x - y)^2
}


############################################################
# Model Confidence Set
############################################################

MCS_f <- function(LOSS, B, alpha) {

  # Computes the Model Confidence Set using a moving-block
  # bootstrap.
  #
  # Inputs:
  #   LOSS  : T x nmod matrix of model-specific loss values,
  #           where T is the number of observations and nmod
  #           is the number of competing models;
  #   B     : number of bootstrap replications;
  #   alpha : significance level of the MCS procedure.
  #
  # Output:
  #   Object returned by MCSprocedure::mcsTest(), containing
  #   the Model Confidence Set results and associated test
  #   statistics.
  #
  # The block length is selected automatically using b.star(),
  # following Politis and White (2004) with the Patton et al.
  # correction, and is constrained to be at least 2.

  k_len <- max(
    b.star(LOSS)[, 2]
  )

  k_len <- ifelse(
    k_len < 1.99,
    2,
    k_len
  )

  temp <- mcsTest(
    LOSS,
    alpha = alpha,
    nboot = B,
    nblock = k_len,
    boot = c("block")
  )

  return(temp)
}