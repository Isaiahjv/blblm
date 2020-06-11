#' @import purrr
#' @import furrr
#' @import stats
#' @import future
#' @import parallel
#' @import MASS
#' @import utils
#' @importFrom magrittr %>%
#' @aliases NULL
#' @export blblm
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))



#' Retruns model fitted with bootstrap sampling
#'
#' @param formula model specification of response and predictor variables.
#'
#' @param data dataset to fit
#' @param m integer
#' @param B integer
#' @param model linear model
#' @param para logical indicating parallel computing
#' @param core integer number of cores for paralle
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000, model = 'lr', para = FALSE, core = 2) {
  #browser()
  if (para == FALSE){
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = .,model ,n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)

  } else {
    suppressWarnings(future::plan(future::multiprocess, workers = core))
    data_list <- split_data(data, m)
    estimates <-  data_list %>% furrr::future_map(function(dat){
      lm_each_subsample(formula = formula, data = dat, model, n = nrow(data), B = B)
    })
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}


# split data into m parts of approximated equal sizes
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


# compute the estimates
lm_each_subsample <- function(formula, data,model,n, B) {
  replicate(B, lm_each_boot(formula, data, n, model), simplify = FALSE)
}


# compute the regression estimates for a blb dataset
lm_each_boot <- function(formula, data, n, model = 'lr') {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  if (model == 'lr'){
    lm1(formula, data, freqs)
  }
  else if (model == 'poisson'){
    lm2(formula, data, freqs)
  }

  else if(model == 'neg'){
    lm3(formula, data, freqs)
  }

  else{
    "Not an option"
  }

}

# estimate the negative bionomial regression estimates based on given the number of repetitions
lm3 <- function(formula,data,freqs){

  environment(formula) <- environment()
  fit <- MASS::glm.nb(formula, data = data,weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))

}

# estimate the poisson regression estimates based on given the number of repetitions
lm2 <- function(formula,data,freqs){

  environment(formula) <- environment()
  fit <- glm(formula,data,family="poisson",freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))

}

# estimate the linear regression estimates based on given the number of repetitions
lm1 <- function(formula, data, freqs) {

  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#'
#' @param fit linear model
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#'
#' @param fit linear model
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' prints output for model
#'
#' @param x linear model
#'
#' @param ... so on
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' returns bootstrapped estimates of sigma (variance) with optional confidence interval. Parallel computing is also optional
#'
#' @param object linear model
#'
#' @param confidence logical indicating confidence interval
#' @param level numeric indicating level of confidence
#' @param para logical indicating parallel computing
#' @param core integer indicating number of cores for parallel computing
#' @param ... so on
#'
#' @method sigma blblm
#' @export
sigma.blblm <- function(object, confidence = FALSE,
                        level = 0.95, para  = FALSE, core = 4 ,...) {

  if(para == FALSE){
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
  }
  else{
    suppressWarnings(future::plan(future::multiprocess, workers = core))
    est <- object$estimates
    sigma <- mean(furrr::future_map_dbl(est, ~ mean(furrr::future_map_dbl(., "sigma"))))
    if (confidence) {
      alpha <- 1 - 0.95
      limits <- est %>%
        pmap_mean(~ quantile(furrr::future_map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2)),core) %>%
        set_names(NULL)
      return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
    } else {
      return(sigma)
    }
  }
}


#' Returns bootstrapped estimates of paramter coefficients. Parallel computing is optional.
#'
#' @param object linear model
#'
#' @param para logical indicating parallel computing
#' @param core integer indiciating number of cores
#' @param ... so on
#' @method coef blblm
#' @export
coef.blblm <- function(object, para = FALSE, core = 2, ...) {
  est <- object$estimates
  if(para == FALSE){
    map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
  }
  else{
    pmap_mean(est, ~ map_cbind(., "coef") %>% rowMeans(),core)
  }

}

#' Returns confidence intervals for estimated parameter coefficients. Parallel computing is optional.
#'
#' @param object linear model
#'
#' @param parm parameter for confidence interval
#' @param level numeric indicating level of confidence
#' @param para logical indicating  parallel computing
#' @param core integer indicating number of cores in parallel computing
#' @param ... so on
#'
#' @method confint blblm
#' @export
confint.blblm <- function(object, parm = NULL, level = 0.95, para = FALSE, core = 2, ...) {

  if(para == FALSE){
    if (is.null(parm)) {
      parm <- attr(terms(object$formula), "term.labels")
    }
    alpha <- 1 - level
    est <- object$estimates
    out <- map_rbind(parm, function(p) {
      map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
    if (is.vector(out)) {
      out <- as.matrix(t(out))
    }
    dimnames(out)[[1]] <- parm
    out
  }

  else{
    suppressWarnings(future::plan(future::multiprocess, workers = core))
    if (is.null(parm)) {
      parm <- attr(terms(object$formula), "term.labels")
    }
    alpha <- 1 - level
    est <- object$estimates
    out <- pmap_rbind(parm, function(p) {
      pmap_mean(est, ~ furrr::future_map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
    })
    if (is.vector(out)) {
      out <- as.matrix(t(out))
    }
    dimnames(out)[[1]] <- parm
    out

  }


}

#' Returns predictions from fitted boostrap model on test data (data unseen by the fitted model). Parallel computing os optional.
#'
#' @param object linear model
#'
#' @param new_data dataset used for prediction
#' @param confidence logical indicating confidence interval usage
#' @param level numeric indicating level of confidence. see previous.
#' @param para logical indicating parallel cimouting usage
#' @param core integer indicating number of core in parallel computing
#' @param ... so on
#'
#' @method predict blblm
#' @export
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, para = FALSE, core = 2, ...) {

  if(para == FALSE){
    est <- object$estimates
    X <- model.matrix(reformulate(attr(terms.formula(object$formula,data = new_data), "term.labels")), new_data)
    if (confidence) {
      map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
                 apply(1, mean_lwr_upr, level = level) %>%
                 t())
    } else {
      map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
    }
  } else{
    suppressWarnings(future::plan(future::multiprocess, workers = core))
    est <- object$estimates
    X <- model.matrix(reformulate(attr(terms.formula(object$formula,data = new_data), "term.labels")), new_data)
    if (confidence) {
      pmap_mean(est, ~ pmap_cbind(., ~ X %*% .$coef,core) %>%
                 apply(1, mean_lwr_upr, level = level) %>%
                 t(),core)
    } else {
      pmap_mean(est, ~ pmap_cbind(., ~ X %*% .$coef, core) %>% rowMeans(),core)
    }

  }

}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}

pmap_mean <- function(.x, .f, core = 2, ...) {
  suppressWarnings(future::plan(future::multiprocess, workers = core))
  (furrr::future_map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

pmap_cbind <- function(.x, .f,core = 2, ...) {

  suppressWarnings(future::plan(future::multiprocess, workers = core))
  furrr::future_map(.x, .f, ...) %>% reduce(cbind)
}

pmap_rbind <- function(.x, .f, core = 2, ...) {
  suppressWarnings(future::plan(future::multiprocess, workers = core))
  furrr::future_map(.x, .f, ...) %>% reduce(rbind)
}
