library(DAAG)
set.seed(123456)

folder_path <- " "	## directory of the folder

quake_data_all <- read.csv(paste0(folder_path, "quake.csv"))

quake_data_all = c(quake_data_all)$x

T_test = 10

y_train = quake_data_all[-tail(1:length(quake_data_all), T_test)]

y_test = tail(quake_data_all, T_test)

T = length(y_train)

B_constant_grid <- seq(0.900, 0.999, by = 0.001)

y_pred_all = matrix(-9999, T_test, length(B_constant_grid))

MSE = c()

for(grid in 1:length(B_constant_grid))
{
## y_t ~ pois(eta_t)
## eta_t|D_(t-1) ~ gamma(alpha_t, beta_t)
## theta_t| D_t ~ [m_t, C_t]

## log(eta_t) = lambda_t = F' theta_t
## theta_t = G theta_t + w_t
## w_t ~ [0, W_t]

n_state = 2

FF = matrix(c(1, 0), nrow = 1)

GG = matrix(c(1,0,1,1), n_state, n_state)

## theta_0 ~ [m0, C0]
m0 = c(0, 0)
C0 = diag(rep(100, 2))

####
B_constant = B_constant_grid[grid]

B = diag(rep(1/B_constant, n_state))

mm = matrix(-99999, n_state, T + 1)
mm[ , 1] = m0

CC <- list()
CC[[1]] <- C0

for(t in 1:T)
	{
	
	## theta_t|D_(t-1) ~ [a_t, R_t]

	a_t = GG %*% mm[ , t]
	R_t = B %*%  GG %*% CC[[t]] %*% t(GG) %*% B
	
	### log(eta_t)|D_(t-1) ~ [f_t, q_t]
	## s_t is the cov(F theta_t, theta_t)
	
	s_t = R_t %*% t(FF)
	
	f_t = c(FF %*% a_t)
	q_t = c(FF %*% s_t)
	
	### prior of eta_t|D_(t-1) ~ gamma(alpha_t, beta_t)
	
	### E[log (eta_t) | D_(t-1)] = digamma(alpha_t) - log (beta_t) = f_t
	###	Var(log (eta_t) | D_(t-1)) = trigamma(alpha_t) = q_t

#if(2<1){	
	fun = function (x) trigamma(x) - q_t
 	alpha_t = uniroot(fun, lower=0.0001, upper=10, extendInt = "yes")$root
	
	beta_t = exp(digamma(alpha_t) - f_t)
#}

#	alpha_t = 1/q_t
#	beta_t = exp(-f_t)/q_t
	
	### posterior of eta_t | D_t ~ gamma(post_alpha_t, post_beta_t)
	
	post_alpha_t = alpha_t + y_train[t]
	post_beta_t = beta_t + 1
	
	### g_t = E[log (eta_t) | D_t] = digamma(post_alpha_t) - log (post_beta_t)
	### p_t = Var(log (eta_t) | D_t) = trigamma(post_alpha_t)
	
	g_t = digamma(post_alpha_t) - log (post_beta_t)
	p_t = trigamma(post_alpha_t)
	
	### theta_t | D_t ~ [m_t, C_t]
	
	mm[ , t+1] = a_t + s_t %*%(g_t - f_t)/q_t
	CC[[t+1]] = R_t - s_t %*% t(s_t) *(1- p_t/q_t)/q_t
	
	}

### Forecasting part

T_test = length(y_test)

a_tk = matrix(-99999, n_state, T_test + 1)

a_tk[ , 1] = a_t 

R_tk = list()
R_tk[[1]] = R_t	

y_pred = c()

for(j in 1:T_test)
	{
	
	a_tk[ , j+1] = GG %*% a_tk[ , j]
	
	R_tk[[j+1]] = B %*%  GG %*% R_tk[[j]] %*% t(GG) %*% B
	
	###
	s_tk = R_tk[[j+1]] %*% t(FF)
	
	f_tk = c(FF %*% a_tk[ , j+1])
	q_tk = c(FF %*% s_tk)

	### prior of eta_tk

	fun_tk = function (x) trigamma(x) - q_tk
 	alpha_tk = uniroot(fun_tk, lower=0.0001, upper=10, extendInt = "yes")$root
#uniroot(fun, c(.Machine$double.eps, 100))$root
	
	beta_tk = exp(digamma(alpha_tk) - f_tk)
	
	##### simulate eta_tk
	
	eta_tk = rgamma(1, alpha_tk, beta_tk)
	
	#### simulate data
	
	y_pred[j] = rpois(1, eta_tk)

	}
	
MSE[grid] = mean((y_test - y_pred)^2)

y_pred_all[ , grid] = y_pred

}	
	
pdf(paste0(folder_path, "MSE_plot.pdf"))

plot(1:length(MSE), MSE, "o")
 
points(which.min(MSE), MSE[which.min(MSE)], pch = 16, cex = 1, col = "blue")
	
dev.off()
	
	
	
	
	
