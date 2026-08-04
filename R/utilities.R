#' @importFrom stats glm kmeans dbeta pbeta
#' @importFrom betareg betareg

ini_bound <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  ind_M <- as.numeric(m_star_vec > 0)
  Y_mod <- glm(yi_vec ~ m_star_vec + ind_M + x_i_vec + x_i_vec * ind_M + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[seq_len(6)] <- Y_mod_ini

  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  km_res <- kmeans(m_nz, k)
  gp_res <- km_res$cluster
  alpha_0_vec <- numeric(k)
  alpha_1_vec <- numeric(k)
  psi_est_vec <- numeric(k)
  phi_est_vec <- numeric(k)

  for (i in seq_len(k)) {
    m_gp <- m_nz[gp_res == i]
    x_gp <- x_i_nz[gp_res == i]

    beta_mod <- try(betareg(m_gp ~ x_gp, link = "logit"), TRUE)
    error_ind <- inherits(beta_mod, "try-error")
    if (error_ind) {
      alpha_0_vec[i] <- mean(m_gp)
      alpha_1_vec[i] <- 0
      phi_est_vec[i] <- NA
    } else {
      beta_gp <- as.numeric(beta_mod$coefficients$mean)
      alpha_0_vec[i] <- beta_gp[1]
      alpha_1_vec[i] <- beta_gp[2]
      phi_est_vec[i] <- as.numeric(beta_mod$coefficients$precision)
    }

    psi_est_vec[i] <- length(m_gp) / length(m_nz)
  }

  alpha_0_sort <- order(alpha_0_vec, decreasing = TRUE)
  alpha_0_sorted <- alpha_0_vec[alpha_0_sort]
  alpha_1_sorted <- alpha_1_vec[alpha_0_sort]
  psi_est_sorted <- psi_est_vec[alpha_0_sort]


  M_mod <- glm((1 - ind_M) ~ x_i_vec, family = "binomial")
  Del_est <- as.numeric(M_mod$coefficients)
  ini_par[7:8] <- Del_est
  phi_est <- mean(phi_est_vec, na.rm = TRUE)
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11:(10 + k)] <- alpha_0_sorted
  ini_par[(11 + k):(10 + 2 * k)] <- alpha_1_sorted
  ini_par[(11 + 2 * k):(9 + 3 * k)] <- psi_est_sorted[seq_len(k - 1)]

  lb_est <- c(rep(-Inf, 6), rep(-10, 2), 5, 0.5, rep(-10, k), rep(-10, k), rep(0.01, k - 1))
  ub_est <- c(rep(Inf, 6), rep(10, 2), Inf, Inf, rep(10, k), rep(10, k), rep(0.99, k - 1))
  return(list(ini_par, lb_est, ub_est))
}

mediation_effect_cal <-
  function(para_vec,
           x_1,
           x_2,
           confound_mat,
           x4_inter = TRUE,
           x5_inter = TRUE) {
    num_conf <- ncol(confound_mat)
    k <- (length(para_vec) - 9 - 3 * num_conf) / 3
    beta_0 <- para_vec[1]
    beta_1 <- para_vec[2]
    beta_2 <- para_vec[3]
    beta_3 <- para_vec[4]
    beta_4 <- para_vec[5]
    beta_5 <- para_vec[6]
    if (!x4_inter) {
      beta_4<-0
    }
    if (!x5_inter) {
      beta_5<-0
    }
    gamma_0 <- para_vec[7]
    gamma_1 <- para_vec[8]
    phi <- para_vec[9]
    delta <- para_vec[10]

    alpha_0_vec <- para_vec[11:(10 + k)]
    alpha_1_vec <- para_vec[(11 + k):(10 + 2 * k)]

    if (k == 1) {
      psi_vec <- 1
      beta_conf <- para_vec[seq(13, 12 + num_conf)]
      alpha_conf <-
        para_vec[seq(13 + num_conf, 12 + 2 * num_conf)]
      gamma_conf <-
        para_vec[seq(13 + 2 * num_conf, 12 + 3 * num_conf)]
    } else {
      psi_temp <- para_vec[(11 + 2 * k):(9 + 3 * k)]
      psi_vec <- c(psi_temp, 1 - sum(psi_temp))
      beta_conf <- para_vec[seq(10 + 3 * k, 9 + 3 * k + num_conf)]
      alpha_conf <-
        para_vec[seq(10 + 3 * k + num_conf, 9 + 3 * k + 2 * num_conf)]
      gamma_conf <-
        para_vec[seq(10 + 3 * k + 2 * num_conf, 9 + 3 * k + 3 * num_conf)]
    }

    conf_mean <- colMeans(confound_mat)

    delta_1_x2 <-
      gamma_0 + gamma_1 * x_2 + sum(gamma_conf * conf_mean)
    mu_x2 <-
      alpha_0_vec + alpha_1_vec * x_2 + sum(alpha_conf * conf_mean)

    delta_1_x1 <-
      gamma_0 + gamma_1 * x_1 + sum(gamma_conf * conf_mean)
    mu_x1 <-
      alpha_0_vec + alpha_1_vec * x_1 + sum(alpha_conf * conf_mean)

    NIE_1 <-
      (beta_1 + beta_5 * x_2) * (expect_M_x(psi_vec, delta_1_x2, mu_x2) - expect_M_x(psi_vec, delta_1_x1, mu_x1))
    NIE_2 <-
      (beta_2 + beta_4 * x_2) * (expit(delta_1_x1) - expit(delta_1_x2))
    NDE <-
      beta_3 * (x_2 - x_1) + beta_4 * (x_2 - x_1) * (1 - expit(delta_1_x1)) + beta_5 * (x_2 - x_1) * expect_M_x(psi_vec, delta_1_x1, mu_x1)
    NIE <- NIE_1 + NIE_2

    return(c(NIE_1, NIE_2, NDE, NIE))
  }

li_total_raw <- function(para_vec, yi_vec, obs_m_vec, xi_vec, li_vec, confound_mat) {
  lk_value_vec <- numeric(length(yi_vec))
  for (i in seq_len(length(yi_vec))) {
    yi <- yi_vec[i]
    obs_m <- obs_m_vec[i]
    xi <- xi_vec[i]
    li <- li_vec[i]
    confound_vec <- confound_mat[i, ]

    if (obs_m > 1e-50) {
      lk_value_temp <- li_1_raw_func(para_vec, yi, obs_m, xi, confound_vec)
      lk_value_vec[i] <- lk_value_temp
    } else if (obs_m < 1e-50 && obs_m >= 0) {
      lk_value_temp <- li_2_raw_func(para_vec, yi, xi, li, confound_vec)
      lk_value_vec[i] <- lk_value_temp
    } else {
      warning("negative m")
    }
  }
  li_total <- -sum(lk_value_vec)
  return(li_total)
}

ini_bound_nomix <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  ind_M <- as.numeric(m_star_vec > 0)
  Y_mod <- glm(yi_vec ~ m_star_vec + ind_M + x_i_vec + x_i_vec * ind_M + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[seq_len(6)] <- Y_mod_ini

  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  beta_mod <- try(betareg(m_nz ~ x_i_nz, link = "logit"), TRUE)
  error_ind <- inherits(beta_mod, "try-error")
  if (error_ind) {
    alpha_0 <- mean(m_nz)
    alpha_1 <- 0
    phi_est <- NA
  } else {
    beta_gp <- as.numeric(beta_mod$coefficients$mean)
    alpha_0 <- beta_gp[1]
    alpha_1 <- beta_gp[2]
    phi_est <- as.numeric(beta_mod$coefficients$precision)
  }
  M_mod <- glm((1 - ind_M) ~ x_i_vec, family = "binomial")
  Del_est <- as.numeric(M_mod$coefficients)

  ini_par[7:8] <- Del_est
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11] <- alpha_0
  ini_par[12] <- alpha_1

  lb_est <- c(rep(-Inf, 6), -10, -10, 0.1, 0.5, -10, -10)
  ub_est <- c(rep(Inf, 6), 10, 10, Inf, Inf, 10, 10)


  return(list(ini_par, lb_est, ub_est))
}

ini_bound_nz <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  Y_mod <- glm(yi_vec ~ m_star_vec + x_i_vec + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[c(1, 2, 4, 6)] <- Y_mod_ini
  ini_par[c(3, 5)] <- 0
  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  km_res <- kmeans(m_nz, k)
  gp_res <- km_res$cluster
  alpha_0_vec <- numeric(k)
  alpha_1_vec <- numeric(k)
  psi_est_vec <- numeric(k)
  phi_est_vec <- numeric(k)

  for (i in seq_len(k)) {
    m_gp <- m_nz[gp_res == i]
    x_gp <- x_i_nz[gp_res == i]

    beta_mod <- try(betareg(m_gp ~ x_gp, link = "logit"), TRUE)
    error_ind <- inherits(beta_mod, "try-error")
    if (error_ind) {
      alpha_0_vec[i] <- mean(m_gp)
      alpha_1_vec[i] <- 0
      phi_est_vec[i] <- NA
    } else {
      beta_gp <- as.numeric(beta_mod$coefficients$mean)
      alpha_0_vec[i] <- beta_gp[1]
      alpha_1_vec[i] <- beta_gp[2]
      phi_est_vec[i] <- as.numeric(beta_mod$coefficients$precision)
    }

    psi_est_vec[i] <- length(m_gp) / length(m_nz)
  }

  alpha_0_sort <- order(alpha_0_vec, decreasing = TRUE)
  alpha_0_sorted <- alpha_0_vec[alpha_0_sort]
  alpha_1_sorted <- alpha_1_vec[alpha_0_sort]
  psi_est_sorted <- psi_est_vec[alpha_0_sort]

  Del_est <- c(-100, 0)
  ini_par[7:8] <- Del_est
  phi_est <- mean(phi_est_vec, na.rm = TRUE)
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11:(10 + k)] <- alpha_0_sorted
  ini_par[(11 + k):(10 + 2 * k)] <- alpha_1_sorted
  ini_par[(11 + 2 * k):(9 + 3 * k)] <- psi_est_sorted[seq_len(k - 1)]

  lb_est <- c(-Inf, -Inf, 0, -Inf, 0, -Inf, -100, 0, 5, 0.5, rep(-10, k), rep(-10, k), rep(0.01, k - 1))
  ub_est <- c(Inf, Inf, 0, Inf, 0, Inf, -100, 0, Inf, Inf, rep(10, k), rep(10, k), rep(0.99, k - 1))

  return(list(ini_par, lb_est, ub_est))
}

ini_bound_nz_nomix <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  Y_mod <- glm(yi_vec ~ m_star_vec + x_i_vec + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[c(1, 2, 4, 6)] <- Y_mod_ini
  ini_par[c(3, 5)] <- 0

  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  beta_mod <- try(betareg(m_nz ~ x_i_nz, link = "logit"), TRUE)
  error_ind <- inherits(beta_mod, "try-error")
  if (error_ind) {
    alpha_0 <- mean(m_nz)
    alpha_1 <- 0
    phi_est <- NA
  } else {
    beta_gp <- as.numeric(beta_mod$coefficients$mean)
    alpha_0 <- beta_gp[1]
    alpha_1 <- beta_gp[2]
    phi_est <- as.numeric(beta_mod$coefficients$precision)
  }
  Del_est <- c(-100, 0)

  ini_par[7:8] <- Del_est
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11] <- alpha_0
  ini_par[12] <- alpha_1

  lb_est <- c(-Inf, -Inf, 0, -Inf, 0, -Inf, -100, 0, 0.1, 0.5, -10, -10)
  ub_est <- c(Inf, Inf, 0, Inf, 0, Inf, -100, 0, Inf, Inf, 10, 10)

  return(list(ini_par, lb_est, ub_est))
}

expect_M_x <- function(psi, delta_1, mu) {
  res <- sum(psi * (1 - expit(delta_1)) * expit(mu))
  return(res)
}

expit <- function(x) {
  exp(x) / (1 + exp(x))
}

li_1_raw_func <- function(para_vec, yi, m_star, x_i, confound_vec) {
  num_conf <- length(confound_vec)
  k <- (length(para_vec) - 9 - 3 * num_conf) / 3

  beta_0 <- para_vec[1]
  beta_1 <- para_vec[2]
  beta_2 <- para_vec[3]
  beta_3 <- para_vec[4]
  beta_4 <- para_vec[5]
  beta_5 <- para_vec[6]

  gamma_0 <- para_vec[7]
  gamma_1 <- para_vec[8]
  phi <- para_vec[9]
  delta <- para_vec[10]

  alpha_0_vec <- para_vec[11:(10 + k)]
  alpha_1_vec <- para_vec[(11 + k):(10 + 2 * k)]

  if (k == 1) {
    psi_vec <- 1
    beta_conf <- para_vec[seq(13, 12 + num_conf)]
    alpha_conf <-
      para_vec[seq(13 + num_conf, 12 + 2 * num_conf)]
    gamma_conf <-
      para_vec[seq(13 + 2 * num_conf, 12 + 3 * num_conf)]
  } else {
    psi_temp <- para_vec[(11 + 2 * k):(9 + 3 * k)]
    psi_vec <- c(psi_temp, 1 - sum(psi_temp))
    beta_conf <- para_vec[seq(10 + 3 * k, 9 + 3 * k + num_conf)]
    alpha_conf <-
      para_vec[seq(10 + 3 * k + num_conf, 9 + 3 * k + 2 * num_conf)]
    gamma_conf <-
      para_vec[seq(10 + 3 * k + 2 * num_conf, 9 + 3 * k + 3 * num_conf)]
  }



  mu_1 <-
    expit(alpha_0_vec + alpha_1_vec * x_i + sum(alpha_conf * confound_vec))
  delta_i <-
    expit(gamma_0 + gamma_1 * x_i + sum(gamma_conf * confound_vec))
  norm_part <-
    (
      yi - beta_0 - beta_1 * m_star - beta_2 - (beta_3 + beta_4) * x_i - beta_5 * x_i * m_star - sum(beta_conf *
                                                                                                       confound_vec)
    )^2 / (2 * delta^2)
  beta_part <-
    log((1 - delta_i) * sum(psi_vec * dbeta(m_star, mu_1 * phi, (1 - mu_1) * phi)))
  lh_value <-
    -0.5 * log(2 * pi) - log(delta) - norm_part + beta_part
  return(lh_value)
}

li_2_raw_func <- function(para_vec, yi, x_i, l_i, confound_vec) {
  num_conf <- length(confound_vec)
  k <- (length(para_vec) - 9 - 3 * num_conf) / 3
  beta_0 <- para_vec[1]
  beta_1 <- para_vec[2]
  beta_2 <- para_vec[3]
  beta_3 <- para_vec[4]
  beta_4 <- para_vec[5]
  beta_5 <- para_vec[6]

  gamma_0 <- para_vec[7]
  gamma_1 <- para_vec[8]
  phi <- para_vec[9]
  delta <- para_vec[10]

  alpha_0_vec <- para_vec[11:(10 + k)]
  alpha_1_vec <- para_vec[(11 + k):(10 + 2 * k)]

  if (k == 1) {
    psi_vec <- 1
    beta_conf <- para_vec[seq(13, 12 + num_conf)]
    alpha_conf <-
      para_vec[seq(13 + num_conf, 12 + 2 * num_conf)]
    gamma_conf <-
      para_vec[seq(13 + 2 * num_conf, 12 + 3 * num_conf)]
  } else {
    psi_temp <- para_vec[(11 + 2 * k):(9 + 3 * k)]
    psi_vec <- c(psi_temp, 1 - sum(psi_temp))
    beta_conf <- para_vec[seq(10 + 3 * k, 9 + 3 * k + num_conf)]
    alpha_conf <-
      para_vec[seq(10 + 3 * k + num_conf, 9 + 3 * k + 2 * num_conf)]
    gamma_conf <-
      para_vec[seq(10 + 3 * k + 2 * num_conf, 9 + 3 * k + 3 * num_conf)]
  }

  mu_2 <-
    expit(alpha_0_vec + alpha_1_vec * x_i + sum(alpha_conf * confound_vec))
  delta_i <-
    expit(gamma_0 + gamma_1 * x_i + sum(gamma_conf * confound_vec))
  norm_part <-
    exp(-(yi - beta_0 - (beta_3) * x_i - sum(beta_conf * confound_vec))^2 / (2 * delta^
                                                                               2))

  int_const_0 <-
    exp(-(
      yi - beta_0 - beta_2 - (beta_3 + beta_4) * x_i - sum(beta_conf * confound_vec)
    )^2 / (2 * delta^2))

  int_const_li <-
    exp(-(
      yi - beta_0 - beta_1 / l_i - beta_2 - (beta_3 + beta_4) * x_i - beta_5 * x_i / l_i - sum(beta_conf *
                                                                                                 confound_vec)
    )^2 / (2 * delta^2))

  const_part <- ((int_const_0 + int_const_li) / 2) * (1 - delta_i)

  int_part <-
    sum(psi_vec * pbeta(1 / l_i, mu_2 * phi, (1 - mu_2) * phi))

  # beta_part<-log((1-delta_i)*sum(psi_vec*dbeta(mu_1*phi,(1-mu_1)*phi)))
  lh_value <-
    -0.5 * log(2 * pi) - log(delta) + log(delta_i * norm_part + const_part * int_part)
  return(lh_value)
}



