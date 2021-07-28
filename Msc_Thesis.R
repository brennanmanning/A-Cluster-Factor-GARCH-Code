## ------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(dtw)
library(jmotif)
library(rugarch)  # GARCH
library(rmgarch)  # DCC
library(mgarchBEKK) # BEKK


## ------------------------------------------------------------------------------------------------------------------------------------------------
df <- read_csv("SP500Data.csv", na = "#N/A", guess_max = 6000, col_types = cols(`1/3/2000` = col_skip()))


## ------------------------------------------------------------------------------------------------------------------------------------------------
Returns <- as_tibble(log(select(df, -DATE)) - log(lag(select(df, -DATE))))
Returns <- Returns %>% add_column(date = df$DATE)
Returns <- Returns[-1,]
df_names <- colnames(Returns)
df_names[-length(df_names)] <- substr(df_names[-length(df_names)],1, nchar(df_names[-length(df_names)])-3)
colnames(Returns) <- df_names


## ------------------------------------------------------------------------------------------------------------------------------------------------
full_stocks <- colSums(is.na(Returns)) == 0
Returns <- Returns[,full_stocks]
dates <- select(Returns, date)
Returns <- select(Returns, -date)
Returns <- Returns * 100


## ------------------------------------------------------------------------------------------------------------------------------------------------
paa_returns <- function(x){
    paa(x, 1768)
}

Returns_paa <- Returns %>% map_dfc(paa_returns)


## ------------------------------------------------------------------------------------------------------------------------------------------------
simil_matrix <- function(ts){
    k <- length(ts) # Number of time series
    D <- matrix(NA, nrow = k, ncol = k)
    for (i in 1:k){
        for (j in i:k){
            print(c(i,j))
            D[i,j] <- dtw(ts[,i], ts[,j], distance.only = T)$normalizedDistance
        }
    }
    return(D)
}

c_Medoids <- function(c, m, D, tol, maxiter){
    k <- dim(D)[1] # Number of time series
    init_idx <- sort(sample(1:k, c))
    # Initialize Membership Matrix
    U_old <- matrix(NA, nrow = c, ncol = k)
    U_new <- calculate_membership(init_idx, m, D)
    medoids <- calculate_medoids(U_new, D)
    iter <- 0
    while(all.equal(U_old, U_new, tolerance = tol) != T && iter < maxiter){
        U_old <- U_new
        U_new <- calculate_membership(medoids, m, D)
        medoids <- calculate_medoids(U_new, D)
        iter <- iter + 1
    }
    
    res <- list(start = init_idx, membership = U_new, medoids = medoids, iterations = iter)
    return(res)
}

calculate_membership <- function(medoids, m, D){
    k <- dim(D)[1]
    c <- length(medoids)
    U <- matrix(NA, nrow = c, ncol = k)
    for (i in 1:c){
        for (j in 1:k){
            med <- medoids[i]
            U[i,j] <- D[j,med]^(-1/(m-1))/(colSums(D[medoids,]^(-1/(m-1)))[j])
        }
    }
    U[, medoids] <- diag(c)
    return(U)
}

calculate_medoids <- function(U, D){
    c <- dim(U)[1]
    J <- rowSums(U[1,] * D)
    medoids <- which.min(J)
    if (c > 1){
        for(i in 2:c){
            J <- rowSums(U[i,] * D)
            idx <- which.min(J[-medoids])
            medoids <- c(medoids, idx + sum(medoids < idx))
        }
    }    
    return(medoids)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
D_mat <- simil_matrix(Returns)


## ------------------------------------------------------------------------------------------------------------------------------------------------
D_paa <- simil_matrix(Returns_paa)

D_paa[lower.tri(D_paa, diag = F)] <- t(D_paa)[lower.tri(D_paa, diag = F)]

## ------------------------------------------------------------------------------------------------------------------------------------------------
vol <- read_csv("Volume.csv", na = "#N/A", guess_max = 6000)
vol <- as_tibble(vol)
Volumes <- vol[,full_stocks]
Volumes <- select(Volumes, -DATE)
sort(colSums(Volumes), decreasing = TRUE)[1:10]

## ------------------------------------------------------------------------------------------------------------------------------------------------
estimate_factor_GARCH <- function(X, factors){
    # Get dimensions of model
    Time <- length(X)
    N <- ifelse(is.matrix(X), dim(X)[2], 1)
    k <- dim(factors)[2]
    # Specify univariate model for each factor 
    f_garch_params <- matrix(0, nrow = k, ncol = 4)
    f_garch_std <- matrix(0, nrow = Time, ncol = k)
    ugarch.spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
                                   variance.model = list(model = 'sGARCH', garchOrder = c(1,1)),
                                   distribution.model = 'norm')
    GARCH_fit <- c()
        
    # Get the parameters for each of the 
    for (i in 1:k){
        ugarch.fit <- ugarchfit(ugarch.spec, factors[,i])
        f_garch_params[i,] <- coef(ugarch.fit)
        f_garch_std[,i] <- ugarch.fit@fit$sigma
        GARCH_fit <- c(GARCH_fit, ugarch.fit)
    }
    
    
    start_pars <- c(rep(1,k),1)
    Lambda <- matrix(0,nrow = N, ncol = k)
    Omega <- matrix(0, nrow = N, ncol = N)
    mess <- c()
    for (i in 1:N){
        res_optim <- optim(par = start_pars,
                           fn = function(pars) - likelihood_factors(pars, X[,i], factors, f_garch_std, f_garch_params)$LN,
                           gr = function(pars) - likelihood_factors(pars, X[,i], factors, f_garch_std, f_garch_params)$GN,
                           method = 'BFGS')
        Lambda[i,] <- res_optim$par[1:k]
        Omega[i,i] <- res_optim$par[k+1]
        mess <- c(mess, res_optim$convergence)
    }
    
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            Omega[i,j] <- mean(X[,i] * X[,j]) - sum(apply(Lambda[c(i,j),], 2, prod) * colMeans(factors^2))
        }
    }
    Omega[lower.tri(Omega, diag = F)] <- t(Omega)[lower.tri(Omega, diag = F)]
    return(list(GARCH_params = f_garch_params, GARCH_fit =GARCH_fit, Sigma = f_garch_std, Lambda = Lambda, Omega = Omega, message = mess))
}


likelihood_factors <- function(pars, x, factors, Sigma, ugarch.params){
    Time <- length(x)
    k <- dim(factors)[2]
    lambda <- pars[1:k]
    omega <- pars[k+1]
    alpha <- ugarch.params[,2]
    beta <- ugarch.params[,3] 
    
    likelihood = 0
    grad_omega = 0
    grad_lambda = rep(0,k)
    for (t in 2:Time){
        
        H = alpha * factors[t-1,]^2 + beta * Sigma[t-1,]^2
        h = omega + sum(lambda^2 * H)
            
        likelihood = likelihood - 1/2 * log(2 * pi) - 1/2 * log(abs(h)) - 1/2 * x[t]^2/h
        grad_omega = grad_omega + (-1/2 + x[t]^2/h)/h
        grad_lambda = grad_lambda + (x[t]^2 / h -1)/h * diag(H)  %*% lambda
    }
    
    
    return(list(LN = likelihood, GN = c(grad_lambda, grad_omega)))
}



## ------------------------------------------------------------------------------------------------------------------------------------------------
IntradayReturns <- read_csv("D:/Thesis/IntradayReturns.csv", 
                            col_types = cols(X1 = col_skip(), Datetime = col_datetime(format = "%Y-%m-%d %H:%M:%S")))


## ------------------------------------------------------------------------------------------------------------------------------------------------
intraday_dates <- unique(date(IntradayReturns$Datetime))
return_dates <- mdy(dates$date)
Returns$Date <- return_dates 
Returns_forecast <- filter(Returns, Date %in% intraday_dates)


## ------------------------------------------------------------------------------------------------------------------------------------------------
# Get 2332 days of forecast and thus use 2559 to estimate parameters
# factors is assumed to have last available data point so you want to 
# forecast n.ahead days from that last observation
factor_GARCH_forecast <- function(n.ahead, uGARCH_params, Sigma, factors, Lambda, Omega){
    alpha_0 <- uGARCH_params[,1]
    alpha_1 <- uGARCH_params[,2]
    beta <- uGARCH_params[,3]
    
    forecast_factors <- alpha_0 + alpha_1 * factors^2 + beta * Sigma
    if (n.ahead > 1){
        for (i in 2:n.ahead){
            forecast_factors <- alpha_0 + (alpha_1 + beta) * forecast_factors 
        }
    }
    Sigma_forecast <- diag(forecast_factors)
    forecast <- Lambda  %*% Sigma_forecast  %*% t(Lambda) + Omega
    return(forecast)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
get_fgarch_forecasts <- function(X, factors, n.ahead){
    cov_forecast <- c()
    for (t in 0:2331){
        if (t %% 22 == 0){
            fGARCH <- estimate_factor_GARCH(X[(t+1):(t+2559),], factors[(t+1):(t+2559),])
            Lambda <- fGARCH$Lambda
            Omega <- fGARCH$Omega
            uGARCH_params <- fGARCH$GARCH_params
            Sigma <- fGARCH$Sigma[t+2559,]^2 
        }
        forecast <- factor_GARCH_forecast(n.ahead, uGARCH_params, Sigma, factors[t+2559,], Lambda, Omega)
        cov_forecast <- c(cov_forecast, as.vector(forecast))
    }
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
estimate_RM <- function(X){
    Time <- dim(X)[1]
    H <- X[1,]  %*% t(X[1,])
    lambda <- 0.96
    for (t in 2:Time){
        H <- (1-lambda) * X[t,]  %*% t(X[t,]) + lambda * H
    }
    return(H)
}

forecast_RM <- function(n.ahead, X, H){
    lambda <- 0.96
    H <- (1 - lambda) * X  %*% t(X) + lambda * H
    return(H)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
forecast_BEKK <- function(X, n.ahead, C, A, G, H){
    
    forecast <- t(C)  %*% C + t(A)  %*% X  %*% t(X)  %*% A + t(G)  %*% H  %*% G
    if (n.ahead > 1){
        for (i in 2:n.ahead){
            forecast <- t(C)  %*% C + t(A + G)  %*% forecast  %*% (A + G)
        }
    }
    return(forecast)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
forecast_GOGARCH <- function(n.ahead, X, A, uGARCH_params, Sigma){
    alpha_0 <- uGARCH_params[1,]
    alpha_1 <- uGARCH_params[2,]
    beta <- uGARCH_params[3,]
    
    forecast_factors <- alpha_0 + alpha_1 * X^2 + beta * Sigma 
    if (n.ahead > 1){
        for (i in 2:n.ahead){
            forecast_factors <- alpha_0 + (alpha_1 + beta) * forecast_factors
        }
    }
    
    forecast <- A  %*% diag(forecast_factors)  %*% t(A)
    return(forecast)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
forecast_DCC <- function(n.ahead, X, dcc_params, Qbar, Q, uGARCH_params, Sigma){
    alpha_0 <- uGARCH_params[1,]
    alpha_1 <- uGARCH_params[2,]
    beta <- uGARCH_params[3,]
    
    D <- alpha_0 + alpha_1 * X^2 + beta * Sigma
    if (n.ahead > 1){
        for ( i in 2:n.ahead){
            D <- alpha_0 + (alpha_1 + beta) * D
        }
    }
    
    D <- diag(D)
    
    Q1 <- (1 - sum(dcc_params)) * Qbar + dcc_params[1] * X  %*% t(X) + dcc_params[2] * Q
    Qk <- (1 - sum(dcc_params)^(n.ahead - 1)) * Qbar + sum(dcc_params)^(n.ahead - 1) * Q1
    
    Qdiag <- sqrt(diag(Qk))
    
    Rk <- solve(diag(Qdiag))  %*% Qk  %*% solve(diag(Qdiag))
    
    H <- D  %*% Rk  %*% D
    return(H)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
get_other_forecasts <- function(X, n.ahead){
    N <- dim(X)[2]
    RM_forecasts <- c()
    DCC_forecasts <- c()
    #GOGARCH_forecasts <- c()

    
    ugarch_spec <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1,1)), 
                              mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                             distribution.model = 'norm')
    uspec <- multispec(replicate(N, ugarch_spec))
    dcc <- dccspec(uspec, model = 'DCC', distribution = 'mvnorm')
    #gogarch <- gogarchspec(mean.model = list(model = 'constant'), variance.model = list(model = 'sGARCH', garchOrder = c(1,1)),
                           #distribution.model = 'mvnorm' )
    lambda <- 0.96
    
    for (t in 0:2331){
        if (t %% 22 == 0){
            dcc_fit <- dccfit(dcc, X[(t+1):(t+2559),])
            #gogarch_fit <- gogarchfit(gogarch, X[(t+1):(t+2559),])
            
            H_last_RM <- estimate_RM(X[(t+1):(t+2559),])
            
            #ugogarch_params <- coef(gogarch_fit)
            #A_gogarch <- as.matrix(gogarch_fit)
            #Sigma_gogarch <- gogarch_fit@mfit$factor.sigmas[2559,]^2
            
            ugarchdcc_params <- as.matrix(coef(dcc_fit, type = 'garch'))
            dim(ugarchdcc_params) <- c(3, N)
            dcc_params <- as.vector(coef(dcc_fit, type = 'dcc'))
            Sigma_dcc <- as.vector(sigma(dcc_fit)[2559,])^2
            Qbar <- dcc_fit@mfit$Qbar
            Q_last <- dcc_fit@mfit$Q[[2559]]
            
            print(t)
        }
        else {
            
            H_last_RM <- (1 - lambda) * X[t+2558,]  %*% t(X[t+2558,]) + lambda * H_last_RM
            
            #Sigma_gogarch <- ugogarch_params[1,] + ugogarch_params[2,] * X[t+2558,]^2 + ugogarch_params[3,] * Sigma_gogarch
            
            Sigma_dcc <- ugarchdcc_params[1,] + ugarchdcc_params[2,] * X[t+2558,]^2 + ugarchdcc_params[3,] * Sigma_dcc
            Q_last <- (1 - sum(dcc_params)) * Qbar + dcc_params[1] * X[t+2558,]  %*% t(X[t+2558,]) + dcc_params[2] * Q_last
        }
        
        RM_forecast <- (1 - lambda) * X[t+2559,]  %*% t(X[t+2559,]) + lambda * H_last_RM
        RM_forecasts <-  cbind(RM_forecasts, as.vector(RM_forecast))
        
        #GOGARCH_forecast <- forecast_GOGARCH(n.ahead, X[t+2559,], A_gogarch, ugogarch_params, Sigma_gogarch)
        #GOGARCH_forecasts <- cbind(GOGARCH_forecasts, as.vector(GOGARCH_forecast))
        
        DCC_forecast <- forecast_DCC(n.ahead, X[t+2559,], dcc_params, Qbar, Q_last, ugarchdcc_params, Sigma_dcc)
        DCC_forecasts <- cbind(DCC_forecasts, as.vector(DCC_forecast))
    }
    
    return(list(RM = RM_forecasts, DCC = DCC_forecasts))
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
demean <- function(x){
    return(x - mean(x))
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
portfolio <- Returns_forecast %>% select(AAPL, AMD, BAC, CSCO, GE, MSFT, MU, T, WFC)
portfolio <- portfolio %>%  mutate_all(demean)


## ------------------------------------------------------------------------------------------------------------------------------------------------
other_forecasts_5days <- get_other_forecasts(as.matrix(portfolio), 5)


## ------------------------------------------------------------------------------------------------------------------------------------------------
IntradayRV <- read_csv("IntradayRV.csv", col_types = cols(X1 = col_skip()))


## ------------------------------------------------------------------------------------------------------------------------------------------------
IntradayRV <- t(IntradayRV)


## ------------------------------------------------------------------------------------------------------------------------------------------------
vech_estimates <- function(cov){
    Time <- dim(cov)[2]
    k <- dim(cov)[1]
    
    vech_mat <- matrix(nrow = Time, ncol = (k + sqrt(k))/2)
    
    for (t in 1:Time){
        X <- matrix(cov[,t], nrow = 9, ncol = 9)
        vech_mat[t,] <- vech(X)
    }
    
    return(vech_mat)
}



get_loss_value(estimate, rv){
    val <- t(estimate - rv)  %*% (estimate - rv)
    return(val)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------
generate_bootstrap_idx <- function(Time, l, B){
    idx <- matrix(nrow = B, ncol = Time)
    n <- ceiling(Time / l)
    
    for (b in 1:B){
        idx_b <- c()
        for(i in 1:n){
            start <- sample(seq(1,Time), size = 1)
            if (start + l - 1 <= Time){
                idx_b <- c(idx_b, seq(start, start + l - 1))
            }
            else {
                runover <- start + l - 1 - Time
                idx_b <- c(idx_b,seq(start, Time))
                idx_b <- c(idx_b, seq(1, runover))
            }
        }
        
        idx[b,] <- idx_b
    }
    
    return(idx)
}

