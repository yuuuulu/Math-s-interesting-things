bs_mixture_em <- function(t, delta, alpha_0, beta_0) {
  # 函数参数说明 ----------------------------------------------------------
  # t: 样本
  # delta: 指示向量 (0=右删失, 1=完全观测)
  # alpha_0: alpha初始值
  # beta_0: beta初始值
  # 返回: 包含最终参数和收敛验证的列表
  
  f_1 <- function(t, alpha, beta) {
    (1/(sqrt(2*pi)*alpha*beta)) * (beta/t)^0.5 * 
      exp(-(1/(2*alpha^2)) * (t/beta + beta/t - 2))
  }
  
  f_2 <- function(t, alpha, beta) {
    (1/(sqrt(2*pi)*alpha*beta)) * (beta/t)^1.5 * 
      exp(-(1/(2*alpha^2)) * (t/beta + beta/t - 2))
  }
  
  calculate_epsilon <- function(t, alpha, beta) {
    numerator <- 0.5 * f_1(t, alpha, beta)
    denominator <- numerator + 0.5 * f_2(t, alpha, beta)
    
  }
  
  ig_cdf <- function(q, alpha, beta) {
    pgamma(beta / q, shape = alpha, scale = 1, lower.tail = TRUE)
  }
  
  rig_cdf <- function(q, alpha, beta) {
    1 - pgamma(beta / q, shape = alpha, scale = 1, lower.tail = TRUE)
  }
  
  S_1 <- function(t, alpha, beta) 1 - ig_cdf(t, alpha, beta)
  S_2 <- function(t, alpha, beta) 1 - rig_cdf(t, alpha, beta)
  
  calculate_gamma <- function(t, alpha, beta) {
    s1 <- S_1(t, alpha, beta)
    s2 <- S_2(t, alpha, beta)
    s1 / (s1 + s2)
  }
  
  
  f_bs_mixture <- function(t, alpha, beta) {
    0.5 * f_1(t, alpha, beta) + 0.5 * f_2(t, alpha, beta)
  }
  
  calculate_T_star <- function(t, alpha, beta) {
    tryCatch({
      #求积分
      integrand <- function(x) x * f_bs_mixture(x, alpha, beta)
      integral <- integrate(integrand, t, Inf, rel.tol = 1e-8)$value
      survival = 1 - pnorm(1/alpha*(sqrt(t/beta)-sqrt(beta/t)))
      ifelse(survival > 1e-10, integral/survival, t)
    }, error = function(e) t)
  }
  
  n <- length(t)
  #10^-5
  tol <- 1e-5
  params_history <- data.frame(
    alpha = numeric(0),
    beta = numeric(0)
  )
  
  # 初始参数
  current <- list(
    alpha = alpha_0,
    beta = beta_0
  )
  params_history <- rbind(params_history, current)
  
  converged <- FALSE
  iter <- 0
  
  while (!converged) {
    iter <- iter + 1
    
    # E步：
    uncensored <- t[delta == 1]
    censored <- t[delta == 0]
    
    epsilon <- sapply(uncensored, calculate_epsilon, 
                      current$alpha, current$beta)
    gamma <- sapply(censored, calculate_gamma, 
                    current$alpha, current$beta)
    T_star <- sapply(censored, calculate_T_star, 
                     current$alpha, current$beta)
    
    # M步：
    sum1 <- sum(gamma) + sum(epsilon)
    sum2 <- sum(T_star) + sum(uncensored)
    sum3 <- sum(1/T_star) + sum(1/uncensored)
    
    # Beta更新
    g <- function(beta) {
      sum1 + (n/2) * (-sum2 + beta^2 * sum3)/(sum2 + beta^2 * sum3) - n/2
    }
    beta_new <- uniroot(g, c(1e-6, 1e6), tol = 1e-12)$root
    
    # Alpha更新
    alpha_sq <- (1/n) * (sum2/beta_new + beta_new*sum3 - 2*n)
    alpha_new <- sqrt(pmax(alpha_sq, 1e-8))
    
    # 记录参数历史
    new_params <- list(alpha = alpha_new, beta = beta_new)
    params_history <- rbind(params_history, new_params)
    
    # 收敛检查（使用最近两次迭代）
    last_two <- tail(params_history, 2)
    alpha_change <- abs(diff(last_two$alpha)) / last_two$alpha[1]
    beta_change <- abs(diff(last_two$beta)) / last_two$beta[1]
    
    if (all(c(alpha_change, beta_change) <= tol)) {
      # 最终验证：检查最后三次迭代的稳定性
      last_three <- tail(params_history, 3)
      final_check <- sapply(last_three, function(x) {
        abs(diff(range(x))) / mean(x)
      })
      
      if (all(final_check <= tol)) {
        converged <- TRUE
      }
    }
    
    current <- new_params
  }
  
  # 最终输出验证
  final_params <- tail(params_history, 1)
  output <- list(
    alpha = final_params$alpha,
    beta = final_params$beta,
    iterations = iter,
    convergence = list(
      criteria = "Relative change < 1e-5",
      last_alpha_change = alpha_change,
      last_beta_change = beta_change,
      params_history = params_history
    )
  )
  
  return(output)
}

# 使用示例（AI给的测试的） ----------------------------------------------------------------
# 生成测试数据
set.seed(123)
true_alpha <- 0.5
true_beta <- 2.0
n_samples <- 100
t <- rweibull(n_samples, shape = 2, scale = true_beta)
delta <- rbinom(n_samples, 1, 0.8)  # 20%的删失率

# 运行算法
result <- bs_mixture_em(t, delta, alpha_0 = 1.0, beta_0 = 1.5)

# 打印结果
cat("最终参数估计:\n",
    "alpha:", result$alpha, "\n",
    "beta:", result$beta, "\n",
    "迭代次数:", result$iterations, "\n")