#### DATA PROCESSING FUNCTIONS
library(rootSolve)
library(purrr)

# aggregate impulse response curve
f <- function(t, t_start = 0, tau = 1/exp(1)){
  ret <- sapply(t, function(x){
    if(x <= t_start) return(0)
    fin <- ((x-t_start)/tau)*exp(1-(x-t_start)/tau)
    return(fin)
  })
  return(ret)
}

# definite integral of response curve, default from 0
ff <- function(a, b, t_start = 0, tau = 1/exp(1)){
  ret <- sapply(seq_along(a), function(x){
    s <- ((max(a[x],t_start)-t_start) + tau)*exp(-(max(a[x],t_start)-t_start)/tau)
    f <- ((max(b[x],t_start)-t_start) + tau)*exp(-(max(b[x],t_start)-t_start)/tau)
    fin <- exp(1)*(s-f)
  })
  return(ret)
}

# sum from all touches of an observation minus the threshold
# used to find zeros of the function
mu_h <- function(t, t_starts, h, tau = 1/exp(1)){
  m <- 0
  n <- length(t_starts)
  for(i in 1:n){
    m <- m + f(t, t_starts[i], tau)
  }
  return(m-h)
}

# find zeros of mu(t)-h with uniroot.all()
t_over_h <- function(mu_h, t_starts, h, tau = 1/exp(1), res_time = 10){
  mu_h2 <- function(x){
    return(mu_h(x, t_starts, h, tau))
  }
  roots <- tryCatch({uniroot.all(mu_h2, c(min(t_starts), res_time))},
                    error = function(e){return(c(0,0))})
  return(roots)
}

# put output of t_over_h into vectors
process_roots <- function(roots, res_time = 10){
  n <- floor((length(roots)+1)/2)
  a_s <- rep(0, n)
  b_s <- a_s
  for(i in 1:(n-1)){
    a_s[i] <- min(roots[2*i-1], res_time)
    b_s[i] <- min(roots[2*i], res_time)
  }
  a_s[n] <- roots[2*n-1]
  if(n == length(roots)/2){
    b_s[n] <- min(roots[2*n], res_time)
  } else {
    b_s[n] <- res_time
  }
  return(list(a = a_s, b = b_s))
}

# aggregated super-threshold impact
# this is M() in the text
int_mu_h <- function(a_s, b_s, t_starts, h, tau = 1/exp(1)){
  if(identical(a_s, 0) & identical(b_s, 0)) return(0)
  int <- 0
  n <- length(t_starts)
  m <- length(b_s)
  for(j in 1:m){
    for(i in 1:n){
      int <- int + ff(a_s[j], b_s[j], t_starts[i])
    }
    int <- int + (a_s[j] - b_s[j])*h
  }
  if(is.na(int)) return(0)
  return(int)
}

# length of time above threshold
ell <- function(a_s, b_s){
  s <- sum(b_s - a_s)
  if(is.na(s)) return(0)
  return(s)
}

# Create our matrix (transpose)
make_M <- function(t_list, h, tau = 1/exp(1), multi = TRUE, ones = TRUE){
  # multi is for multiple channels
  # ones = TRUE if we want intercept in matrix
  if(!multi){
    M <- matrix(1, nrow = 4, ncol = length(t_list[[1]]))
    
    M[2:4,]  <- matrix(unlist(lapply(t_list, function(x){
      if(is.null(x)) return(list(0, 0, 0))
      x <- na.omit(x)
      res <- max(x)
      at_res <- mu_h(res, x, h, tau)
      x <- x[x != res]
      roots <- process_roots(t_over_h(mu_h, x, h, tau, res))
      int <- int_mu_h(roots$a, roots$b, x, h, tau)
      l <- tryCatch({ell(roots$a, roots$b)}, error = function(e)return(0))
      return(list(int = int, at_res = at_res, l = l))
    })), nrow = 3, byrow = FALSE)
    
  } else {
    channels <- names(t_list)
    nrows <- 3*length(channels) + 1
    M <- matrix(1, nrow = nrows, ncol = length(t_list[[1]]))
    for(i in seq_along(channels)){
      M[(3*i-1):(3*i+1),] <- matrix(unlist(lapply(t_list[[channels[i]]], function(x){
        if(is.null(x)) return(list(0, 0, 0))
        x <- na.omit(x)
        res <- max(x)
        at_res <- mu_h(res, x, h, tau)
        x <- x[x != res]
        roots <- process_roots(t_over_h(mu_h, x, h, tau, res))
        int <- int_mu_h(roots$a, roots$b, x, h, tau)
        l <- tryCatch({ell(roots$a, roots$b)}, error = function(e)return(0))
        return(list(int = int, at_res = at_res, l = l))
      })), nrow = 3, byrow = FALSE)
    }
  }
  
  M <- data.frame(t(M))
  names(M) <- c('one', 'int_mu_h', 'f_val', 'ell')
  if(ones) return(M)
  return(M[,2:4])
  
}


### ALGORITHMIC FUNCTIONS

sigmoid <- function(x){
  return(1/(1+exp(-x)))
}

# b are results of each obs i
# M has n rows, 1 for each obs
loss <- function(beta, b, M){
  l <- 0
  n <- length(b)
  y <- sigmoid(beta %*% t(M))
  for(i in 1:n){
    if(b[i]){
      l <- l - log(y[i])
    } else {
      l <- l - log(1-y[i])
    }
  }
  return(l/n)
}

# gradient of loss with respect to beta
grad_beta <- function(beta, b, M){
  n <- length(b)
  y <- sigmoid(beta %*% t(M))
  g <- (1/n)*(t(M) %*% c(y-b))
  return(g)
}

# gradient of loss with respect to threshold
grad_h <- function(b, beta, M, epsilon = 10e-6, h, t_list){
  M <- t(M)
  # get numerical gradients for value of M
  dM_dh <- lapply(seq_along(t_list), function(x){
    # get roots for each from which to calculate numerical derivative
    eps <- process_roots(t_over_h(mu_h, t_list[[x]], h + epsilon)) 
    # get the actual value
    eps_int <- int_mu_h(eps$a, eps$b, t_list[[x]], h)
    # numerical value of derivative of agg impact to h
    d <- (eps_int - M[2,x]) / epsilon # two is location of m value
    # numerical value of derivate of length of l to h
    d_2 <- (ell(eps$a, eps$b) - M[4,x])
    return(c(d, d_2)) # other values zero
  })
  
  trace_grad <- lapply(seq_along(dM_dh), function(x){
    M <- t(M)
    bTm <- sum(beta*M[,x])
    y <- sigmoid(bTm)
    grad_obs <- beta*(y-b[x])*bTm*(1-bTm)/(y*(1-y)) # beta * scalar
    grad_obs <- grad_obs[c(2,4)] # other entires of dM_dh are 0
    return(sum(grad_obs*dM_dh[[x]]))
  })
  
  return(sum(unlist(trace_grad)))
}

# This is the cost that we want to minimize actually
# Runtime is too high due to matrix reconfiguration
cost_h <- function(beta_h, b, t_list, tau = 1/exp(1)){
  beta <- beta_h[1:(length(beta_h)-1)]
  h <- beta_h[length(beta_h)]
  M <- make_M(t_list, h, tau)
  J <- loss(b, beta, M)
  return(J)
}

# the fucnction we use to run the analysis for a set of h values
mca_thresh_optim <- function(b, t_list = NULL, h_s = NULL, tau = 1/exp(1), M_s = NULL){
  # M_s can be pre-processed and in parallel in a list with h
  if(is.null(M_s)){
    if(is.null(h_s)) stop('No threshold values specified.')
    if(is.null(t_list)) stop('No list of timestamps provided.')
    M_s <- lapply(h_s, function(x){
      return(make_M(t_list, x))
    })
  }
  
  n <- ncol(M_s[[1]])
  
  models <- lapply(M_s, function(x){
    beta_guess <- optim(rep(0, n), loss, b = b, M = x)
    return(list(beta = beta_guess$par, exp_beta = exp(beta_guess$par),
                pred = sigmoid(beta_guess$par %*% t(x)), error = beta_guess$value))
  })
  
  return(models)
}

log_reg <- function(b, M, epsilon_grad = 10e-6, 
                    epsilon_stop = 10e-6, eta_beta = .1, 
                    eta_h = 10e-4, max_iter = 10e3){
  
  beta_now <- rep(0, ncol(M))# runif(ncol(M), min = -1, max = 1)
  loss_now <- loss(beta_now, b, M)
  
  i <- 0
  
  while(i < max_iter){
    
    beta_last <- beta_now
    
    beta_now <- c(beta_now - eta_beta*grad_beta(beta_now, b, M))
    
    if(all(abs(beta_now - beta_last) < epsilon_stop)) break
    if(all(abs(beta_now - beta_last) < 10*epsilon_stop)) eta_beta <- eta_beta/2
    
    i <- i + 1
    
  }
  
  return(beta_now)
}

# This is the manual gradient descent
mca_thresh_adam <- function(b, t_list = NULL, h_s = NULL, tau = 1/exp(1), M_s = NULL){
  if(is.null(M_s)){
    if(is.null(h_s)) stop('No threshold values specified.')
    if(is.null(t_list)) stop('No list of timestamps provided.')
    M_s <- lapply(h_s, function(x){
      return(make_M(t_list, x))
    })
  }
  
  n <- ncol(M_s[[1]])
  
  models <- lapply(M_s, function(x){
    beta_guess <- log_reg(b, x)
    return(list(beta = beta_guess, exp_beta = exp(beta_guess),
                pred = sigmoid(beta_guess %*% t(x)), error = loss(beta_guess, b, x)))
  })
  
  return(models)
}

# Plots errors versus h values, output from models in functions above
grid_model_plot <- function(models, h_s){
  errors <- unlist(lapply(models, function(x){
    return(x$error)
  }))
  plot(unlist(h_s), errors)
}

## We need a function that also searches over h's for the correct one

# identifies spieks in a set of signals
id_spikes <- function(signals){
  
  assertthat::assert_that(all(sapply(c('unique_id', 'channel', 'time', 'voltage'),
                                     
                                     function(x) x %in% names(signals))))
  
  keys <- unique(signals$unique_id)
  output <- list()
  
  for(i in keys){
    
    output[[i]] <- list()
    sbset <- signals[signals$unique_id == i,]
    
    for(j in unique(sbset$channel)){
      
      chset <- sbset[sbset$channel == j,]
      assertthat::are_equal(seq(0, nrow(sbset)-1, 1), chset$time)
      
      n <- nrow(chset)
      spikes <- list()
      c <- 1
      
      this_signal <- sbset[sbset$channel == j,]$voltage
      
      for(k in 2:(n-1)){
        if(this_signal[k] > this_signal[k-1] && 
           this_signal[k] > this_signal[k+1]){
          spikes[[c]] <- k
          c <- c + 1
        }
      }
      
      output[[i]][[j]] <- unname(unlist(spikes)) #list of times with activity
    }
  }
  return(output)
}


  