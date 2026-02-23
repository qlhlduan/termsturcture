library(dplyr)
library(ggplot2)
library(quadprog)

compute_grouped_mean <- function(data) {
  n_groups <- nrow(data) %/% 2000
  data %>%
    mutate(group = rep(1:n_groups, each = 2000)) %>%
    group_by(group) %>%
    summarise(across(everything(), mean)) %>%
    select(-group)
}

compute_grouped_cov <- function(data, group_size = 2000) {
  n_groups <- nrow(data) %/% group_size
  data %>%
    mutate(
      group = rep(1:n_groups, each = group_size, length.out = n())
    ) %>%
    group_by(group) %>%
    do({
      cov_mat <- cov(select(., -group)) 
      as.data.frame(cov_mat)
    }) %>%
    ungroup() %>%
    select(-group)
}

compute_grouped_sd <- function(data) {
  n_groups <- nrow(data) %/% 2000
  data %>%
    mutate(group = rep(1:n_groups, each = 2000)) %>%
    group_by(group) %>%
    summarise(across(everything(), sd)) %>%
    select(-group)
}


yield <- read.csv("yield.csv")
data <- yield[yield$Date >= "1985-01-01" & yield$Date <= "2000-12-31",]
data <- data[121:181,c(3,6,12,14,19)]
rownames(data) <- 1:nrow(data)

res_bdns <- read.csv("BDNS.csv")[,1:5]
mean_ns <- compute_grouped_mean(res_bdns)
cov_ns <- compute_grouped_cov(res_bdns)

res_aniso <- read.csv("Aniso.csv")[,1:5]
mean_aniso <- compute_grouped_mean(res_aniso)
cov_aniso <- compute_grouped_cov(res_aniso)

res_st <- read.csv("Spatemp.csv")[,1:5]
mean_st <- compute_grouped_mean(res_st)
cov_st <- compute_grouped_cov(res_st)

res_sf <- read.csv("Stat.csv")[,1:5]
mean_sf <- compute_grouped_mean(res_sf)
cov_sf <- compute_grouped_cov(res_sf)

res_nsf <- read.csv("Nonstat.csv")[,1:5]
mean_nsf <- compute_grouped_mean(res_nsf)
cov_nsf <- compute_grouped_cov(res_nsf)

res_gamma <- read.csv("Gamma.csv")[,1:5]
mean_gamma <- compute_grouped_mean(res_gamma)
cov_gamma <- compute_grouped_cov(res_gamma)

res_log <- read.csv("Lognormal.csv")[,1:5]
mean_log <- compute_grouped_mean(res_log)
cov_log <- compute_grouped_cov(res_log)

generate_pmat <- function(cov_matrix, mean_matrix, data, zeta = 4) {
  pmat <- matrix(0, nrow = 61, ncol = 4)
  for (i in 1:61) {
    pvec <- numeric(4)
    for (j in 1:4) {
      # Construct covariance matrix components
      Dmat <- zeta * matrix(as.numeric(c(
        cov_matrix[(5*i-5+j), j], 
        cov_matrix[(5*i-5+j), 5],
        cov_matrix[(5*i), j], 
        cov_matrix[(5*i), 5]
      )), ncol = 2)
      
      dvec <- as.numeric(c(mean_matrix[i, j], mean_matrix[i, 5]))
      rvec <- as.numeric(c(data[i, j], data[i, 5]))
      Amat <- t(matrix(c(1, 0, 1, -1), ncol = 2))
      bvec <- c(1, -50)
      
      qp_j <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
      pvec[j] <- sum(qp_j$solution * rvec)
    }
    pmat[i, ] <- pvec
  }
  return(pmat)
}

calculate_performance_fees <- function(strategy_pmat, benchmark_pmat, delta) {
  calculate_performance_fee <- function(returns_rsdns, returns_benchmark, delta) {
    returns_rsdns <- returns_rsdns / 100 + 1
    returns_benchmark <- returns_benchmark / 100 + 1
    utility_bench <- sum(returns_benchmark - (delta / (2 * (1 + delta))) * returns_benchmark^2)
    
    F <- uniroot(
      function(F) sum((returns_rsdns - F) - (delta / (2 * (1 + delta))) * (returns_rsdns - F)^2) - utility_bench,
      interval = c(-1, 1)
    )$root
    return(F)
  }
  
  sapply(1:ncol(strategy_pmat), function(k) {
    calculate_performance_fee(strategy_pmat[, k], benchmark_pmat[, k], delta)
  })
}

zeta = 1

pmat_ns <- generate_pmat(cov_ns, mean_ns, data, zeta)
pmat_aniso <- generate_pmat(cov_aniso, mean_aniso, data, zeta)
pmat_st <- generate_pmat(cov_st, mean_st, data, zeta)
pmat_sf <- generate_pmat(cov_sf, mean_sf, data, zeta)
pmat_nsf <- generate_pmat(cov_nsf, mean_nsf, data, zeta)
pmat_gamma <- generate_pmat(cov_gamma, mean_gamma, data, zeta)
pmat_log <- generate_pmat(cov_log, mean_log, data, zeta)

colMeans(pmat_ns)

performance_fee_st <- calculate_performance_fees(pmat_st, pmat_ns, delta = 1)
performance_fee_aniso <- calculate_performance_fees(pmat_aniso, pmat_ns, delta = 1)
performance_fee_sf <- calculate_performance_fees(pmat_sf, pmat_ns, delta = 1)
performance_fee_nsf <- calculate_performance_fees(pmat_nsf, pmat_ns, delta = 1)
performance_fee_gamma <- calculate_performance_fees(pmat_gamma, pmat_ns, delta = 1)
performance_fee_log <- calculate_performance_fees(pmat_log, pmat_ns, delta = 1)

percent_st <- performance_fee_st * 10000 / colMeans(pmat_ns)
percent_aniso <- performance_fee_aniso * 10000 / colMeans(pmat_ns)
percent_sf <- performance_fee_sf * 10000 / colMeans(pmat_ns)
percent_nsf <- performance_fee_nsf * 10000 / colMeans(pmat_ns)
percent_gamma <- performance_fee_gamma * 10000 / colMeans(pmat_ns)
percent_log <- performance_fee_log * 10000 / colMeans(pmat_ns)

performance_fee <- rbind(performance_fee_st, performance_fee_aniso, performance_fee_sf, performance_fee_nsf, performance_fee_gamma, performance_fee_log)
rownames(performance_fee) <- c("st", "aniso", "stat", "nonstat", "gamma", "lognormal")

percent_fee <- rbind(percent_st, percent_aniso, percent_sf, percent_nsf, percent_gamma, percent_log)
rownames(percent_fee) <- c("st", "aniso", "stat", "nonstat", "gamma", "lognormal")

print(percent_fee)
