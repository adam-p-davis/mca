setwd('/Users/adamdavis/Documents/Senior/Thesis/Code/')
source('mca_functions.R')
library(ChannelAttribution)
library(reshape2)
library(purrr)
library(ggplot2)
library(gridExtra)
library(markovchain)
library(stargazer)

n_obs <- 100
n_channels <- 3
n_times <- 10 # no days
ones <- TRUE

set.seed(20)
h_real <- 1
beta_real <- c(-1, .2, .2, .2)
poisson_freq <- runif(n_obs)
dts <- rep(1, n_obs)

simulate <- function(n_obs = 100, n_channels = 3, n_times = 10, ones = TRUE,
                     h_real = 1, beta_real = NULL){
  if(is.null(beta_real)) beta_real = c(-1, runif(n_channels))
  poisson_freq <- runif(n_obs)
  dts <- rep(1, n_obs)
  
  simul_spikes <- list()
  
  for(i in 1:n_obs){
    simul_spikes[[as.character(i)]] <- list()
    for(j in 1:n_channels){
      trts <- which(runif(n_times) < poisson_freq[i]*dts[i])
      simul_spikes[[as.character(i)]][[as.character(j)]] <- trts
    }
  }
  
  M_real <- make_M(simul_spikes %>% transpose(), h_real, ones = ones)
  
  b <- sapply(sigmoid(beta_real %*% as.matrix(t(M_real)) + rnorm(n_obs, sd = 1)), function(x)rbinom(1, 1, x))
  
  return(list(spikes = simul_spikes, b = b))
}

pred_and_cv_rates <- function(simul_spikes, b, ones = TRUE, graph = FALSE, give_beta = FALSE){
  n_obs <- length(simul_spikes)
  system.time(sim <- opt_bh(b, simul_spikes %>% transpose(), 0, 2, h_len_max = 25, grid_h = TRUE, ones = ones))
  model_errors <- unlist(lapply(sim$models, function(x)return(x$error)))
  best_ind <- which.min(model_errors)
  model_best <- sim$models[[best_ind]]
  h_best <- sim$h_s[[best_ind]]
  pred_rate <- sum((model_best$pred > .5) == b)/n_obs
  
  folds <- 10
  n <- length(b)
  cv <- rep(0, n)
  guess <- cv
  order <- sample(1:n, n)
  for(i in 1:folds){
    removed <- order[(n*(i-1)/folds):(n*i/folds)]
    new_model <- unlist(mca_thresh_optim(b[-removed], M_s = list(model_best$M[-removed,])),
                        recursive = FALSE)
    guess[removed] <- sigmoid(rowSums(new_model$beta * model_best$M[removed,]))
    cv[removed] <- b[removed] == (guess[removed] > .5)
  }
  
  cv_rate <- sum(cv)/length(cv)
  
  if(graph){
    dfcv <- data.frame(x = guess, b = as.factor(b))
    cv_plot <- ggplot(dfcv, aes(x, fill = b)) + geom_histogram(bins = 50, alpha = .5) +
      scale_fill_manual(values = c('red', 'blue')) + 
      guides(fill = guide_legend(override.aes = list(alpha = .5))) + 
      scale_alpha_manual(values = c(.5), guide = FALSE) +
      labs(x = 'Predicted value', y = 'Frequency', title = 'Cross-validated predictions')
    # cv_plot
    grid.arrange(grid_model_plot(sim$models, sim$h_s), cv_plot, ncol = 1)
  }
  if(give_beta) return(list(pred = pred_rate, cv = cv_rate, beta = model_best$beta))
  return(list(pred = pred_rate, cv = cv_rate))
  
}


# b <- sigmoid(beta_real %*% as.matrix(t(M_real)) + rnorm(n_obs)) > .5


# system.time(sim <- opt_bh(b, simul_spikes %>% transpose(), 0, 2, h_len_max = 11))
# grid_model_plot(sim$models, sim$h_s)


# grid_model_plot(sim$models, sim$h_s)



# model_best <- sim$models[[13]]



# guess <- cv
# for(i in 1:n){
#   M_h2 <- model_best$M[-i,]
#   new_model <- unlist(mca_thresh_optim(b[-i], M_s = list(M_h2)), recursive = FALSE)
#   guess[i] <- sigmoid(sum(new_model$beta*model_best$M[i,]))
#   cv[i] <- b[i] == (guess[i] > .5)
# }


# M_real
# model_best$M

# we need to convert our data to path data
markov_pred_rate <- function(simul_spikes, b, orders = c(1, 2, 3)){
  n_obs <- length(simul_spikes)
  trt_paths <- data.frame(path = rep(NA, n_obs),
                          conversion = rep(0, n_obs),
                          null = rep(0, n_obs))
  for(i in 1:n_obs){
    molt <- melt(simul_spikes[[as.character(i)]])
    path <- paste(paste(molt$L1[order(molt$value)], ' > '), collapse = ' ')
    trt_paths$path[i] <- trimws(substr(path,1,nchar(path) - 4))
    if(b[i]){
      trt_paths$conversion[i] <- 1
    } else {
      trt_paths$null[i] <- 1
    }
  }
  
  new_paths <- aggregate(cbind(conversion, null) ~ path, trt_paths, FUN = sum)
  new_paths$path <- gsub('  ', ' ', 
                         gsub('1', 'one', 
                              gsub('2', 'two',
                                   gsub('3', 'thr', new_paths$path))))
  
  ordered_paths <- gsub(' > ', ',', 
                        gsub('  ', ' ', 
                             gsub('1', 'one', 
                                  gsub('2', 'two',
                                       gsub('3', 'thr', trt_paths$path)))))
  
  rates_to_return <- list()
  
  for(k in orders){
    # train the model
    M <- markov_model(new_paths, 'path', 'conversion', 
                      var_null = 'null', out_more = TRUE, order = k)
    
    # build square transition matrix for markovchain package
    transition <- acast(M$transition_matrix, channel_from ~ channel_to, 
                        value.var = 'transition_probability')
    transition[is.na(transition)] <- 0
    s <- which(rownames(transition) == '(start)')
    transition <- rbind(transition[1:s,], 
                        rep(0, ncol(transition)), 
                        rep(0, ncol(transition)),
                        transition[(s+1):nrow(transition),])
    transition <- cbind(transition[,1:(s-1)],
                        rep(0, nrow(transition)),
                        transition[,s:ncol(transition)])
    
    transition[s+1,s+1] <- 1
    transition[s+2,s+2] <- 1
    
    rownames(transition)[(s+1):(s+2)] <- c('(conversion)', '(null)')
    colnames(transition)[s] <- '(start)'
    
    statez <- colnames(transition)
    chain <- new('markovchain', states = statez, transitionMatrix = transition, name = paste0('chain', k))
    
    n_simul <- 100
    n_states <- length(statez)
    prob_c <- rep(0, n_states)
    for(i in 1:n_states){
      for(j in 1:n_simul){
        simul <- rmarkovchain(n = 100, object = chain, t0 = statez[i])
        if(simul[length(simul)] == '(conversion)'){
          prob_c[i] <- prob_c[i] + 1
        }
      }
      prob_c[i] <- prob_c[i] / n_simul
    }
    
    names(prob_c) <- statez
    
    ordered_ends <- unname(sapply(ordered_paths, function(x)substr(x, max(1, nchar(x) - (4*k -2)), nchar(x))))
    
    rates_to_return[[k]] <- sum(b == unname(sapply(ordered_ends, function(x)prob_c[which(names(prob_c) == x)] > .5)))/n_obs
    
  }
  return(rates_to_return)
}

m_rates <- markov_pred_rate(simul_spikes)
  
  # system.time(M <- markov_model(new_paths, 'path', 'conversion', 
  #                             var_null = 'null', out_more = TRUE, order = 3))
  # c_path <- gsub(' > ', ',', new_paths$path)
  # ends <- unname(sapply(c_path, function(x)substr(x, max(1, nchar(x) - 10), nchar(x))))
  # to_conv <- M$transition_matrix[M$transition_matrix$channel_to == '(conversion)',]
  # ends2 <- ends[ends %in% to_conv$channel_from]
  # c_probs <- unlist(lapply(ends, function(x) to_conv$transition_probability[to_conv == x]))
  # c_probs > .5
  
  # for_pred <- data.frame(ends = ends, c_probs = c_probs, conv = c_probs >= .5)
  # for_pred
  

  # plot(chain, package = 'diagram', box.size = 0.01)
  


## TESTING AND MAKING PLOTS

sims <- simulate()
pred_and_cv_rates(sims$spikes, sims$b)
markov_pred_rate(sims$spikes, b)

out <- 20
h_s <- seq(0, 2, length.out = out)
errors2 <- data.frame(h = h_s, 
                     pred = rep(0, out),
                     cv = rep(0, out),
                     first_order = rep(0, out),
                     second_order = rep(0, out),
                     third_order = rep(0, out))

for(i in 1:length(h_s)){
  sims <- simulate(h_real = h_s[i], n_obs = 500)
  pcv <- pred_and_cv_rates(sims$spikes, sims$b)
  mpr <- markov_pred_rate(sims$spikes, sims$b)
  errors2[i,2:6] <- c(unlist(pcv), unlist(mpr))
}

rate_melt2 <- melt(errors2, id.vars = 'h')

h_plot2 <- ggplot(rate_melt2, aes(h, value, col = variable)) + geom_line() +
  labs(title = 'Predictive accuracy vs. threshold values',
       y = 'Prediction rate') + 
  scale_color_discrete(name = 'Method', labels = c('Impulse-Response Full Fit',
                                                   'IR Cross-validated',
                                                   'First-order Markov',
                                                   'Second-order M',
                                                   'Third-order M'))


h_plot2

errors

rate_melt <- melt(errors, id.vars = 'h')
h_plot <- ggplot(rate_melt, aes(h, value, col = variable)) + geom_line() +
  labs(title = 'Predictive accuracy vs. threshold values',
        y = 'Prediction rate') + 
  scale_color_discrete(name = 'Method', labels = c('Impulse-Response Full Fit',
                                  'IR Cross-validated',
                                  'First-order Markov',
                                  'Second-order M',
                                  'Third-order M'))


h_plot

out <- 20
h_s <- seq(0, 2, length.out = out)
errors <- data.frame(h = h_s, 
                     pred = rep(0, out),
                     cv = rep(0, out),
                     first_order = rep(0, out),
                     second_order = rep(0, out),
                     third_order = rep(0, out))

for(i in 1:length(h_s)){
  sims <- simulate(h_real = h_s[i], n_obs = 500)
  pcv <- pred_and_cv_rates(sims$spikes, sims$b)
  mpr <- markov_pred_rate(sims$spikes, sims$b)
  errors[i,2:6] <- c(unlist(pcv), unlist(mpr))
}

set.seed(50)
h_real <- runif(1)
reps <- 100
beta_real <- c(-1, runif(3))


betas <- matrix(rep(0, reps*4), nrow = reps)
for(i in 1:reps){
  sims <- simulate(beta_real = beta_real, h_real = h_real)
  model <- unlist(mca_thresh_optim(sims$b, sims$spikes %>% transpose(), h_real, ones = TRUE), recursive = FALSE)
  betas[i,] <- model$beta
}

betas_df <- as.data.frame(betas)
names(betas_df) <- c('beta_0', 'beta_1', 'beta_2', 'beta_3')
stargazer(betas_df)

beta_real
h_real
