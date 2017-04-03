library(eegkitdata)
library(ggplot2)
library(gridExtra)
library(stargazer)

setwd('/Users/adamdavis/Documents/Senior/Thesis/Code/')
source('mca_functions.R')

full_wd <- '/Users/adamdavis/Documents/Senior/Thesis/Code/EEGData/eeg_full/'

# eeg_full <- geteegdata(full_wd, cond = 'S1', filename = 'fullS1')

# data(eegdata)
# eegdata <- within(eegdata, {
#   unique_id <- paste(subject, trial, sep = '_')
#   b <- group == 'a'
# })

# load in data that was processed in parallel with poor 
load('EEGData/eeg_rda/eeg_full_S1.rda')
full_1 <- eegdata
load('EEGData/eeg_rda/eeg_full4_S1.rda')
full_2 <- eegdata
load('EEGData/eeg_rda/fullS1 copy.rda')
full_3 <- eegdata
load('EEGData/eeg_rda/fullS1.rda')
full_4 <- eegdata


eegdata <- rbind(full_1, full_2, full_3, full_4)
length(unique(eegdata$subject)) == 122 # number of subjects

# take average of the signals
eeg_agg <- aggregate(voltage ~ subject + channel + time, data = eegdata, FUN = mean)

# identify the points with spikes
spikes <- id_spikes(eeg_agg)

# create output vector
eegdata <- within(eegdata, {
  b <- group == 'a'
})
b <- aggregate(b ~ subject, data = eegdata, FUN = mean)
b2 <- b
b <- b[,2]

# we never want to have to wait that long again
save(eegdata, file = 'full_eegdata.rda')
save(eeg_agg, file = 'full_eegdata_agg.rda')
save(spikes, file = 'spikes.rda')

load('full_eegdata.rda')
load('full_eegdata_agg.rda')
load('spikes.rda')

# this organizes by channel first, how the functions work
# this tranposes the list
spikes <- spikes %>% transpose()

# small testing so we do not waste runtime
smol_spikes <- spikes[1:10]
system.time(res <- opt_bh(b, smol_spikes, 0, 1.5))

# now with large set
# epsilon = default
system.time(results <- opt_bh(b, spikes, 0, 2, h_len_max = 11)) # 23 mins
# epsilon = 1/(2^(h_len+3))
system.time(results2 <- opt_bh(b, spikes, 0, 2, h_len_max = 11)) # 40 mins

# these estimates led us to search over smaller intervals to find the best threshold
system.time(results3 <- opt_bh(b, spikes, .8, 1, h_len_max = 7))

system.time(results4 <- opt_bh(b, spikes, .85, .9, h_len_max = 5))

# put all of these results into 1 to visualize results
size <- 28 # 10 + 10 + 6 + 2
errors <- list(size)
h_s <- errors
ind <- 1
for(i in seq_along(results$models)){
  errors[[ind]] <- results$models[[i]]
  h_s[ind] <- results$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results2$models)){
  errors[[ind]] <- results2$models[[i]]
  h_s[ind] <- results2$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results3$models)){
  errors[[ind]] <- results3$models[[i]]
  h_s[ind] <- results3$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results3$models)){
  errors[[ind]] <- results3$models[[i]]
  h_s[ind] <- results3$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results4$models)){
  errors[[ind]] <- results4$models[[i]]
  h_s[ind] <- results4$h_s[[i]]
  ind <- ind + 1
}

h_s

grid_model_plot(errors, h_s) + 
  labs(title = 'Summary grid of h values from methods with 3 sets of initial parameters') +
  annotate('rect', xmin = .82, xmax = .88, ymin = .36, ymax = .44, alpha = .2, fill = 'blue') +
  annotate('segment', x = .5, xend = .875, y = .5, yend = .4, colour = 'blue', alpha = .5)

# these are the results that we will use
h_best <- .875
model_best <- errors[[4]]
pred_rate <- sum((model_best$pred > .5) == b)/length(model_best$pred)

# leave one out cross validation, we will use the estimated threshold
n <- length(b)
cv <- rep(0, n)
guess <- cv
for(i in 1:n){
  M_h <- model_best$M[-i,]
  new_model <- unlist(mca_thresh_optim(b[-i], M_s = list(M_h)), recursive = FALSE)
  guess[i] <- sigmoid(sum(new_model$beta*model_best$M[i,]))
  cv[i] <- b[i] == (guess[i] > .5)
  
}

n <- length(b)
cv2 <- rep(0, n)
guess2 <- cv
for(i in 1:n){
  M_h2 <- model_best$M[-i,]
  new_model2 <- unlist(mca_thresh_adam(b[-i], M_s = list(M_h2)), recursive = FALSE)
  guess2[i] <- sigmoid(sum(new_model2$beta*model_best$M[i,]))
  cv2[i] <- b[i] == (guess2[i] > .5)
  
}

(cs2_rate <- sum(cv2)/n)


(cv_rate <- sum(cv)/n)
dfcv <- data.frame(x = guess, b = as.factor(b), x2 = guess2)
cv_plot <- ggplot(dfcv, aes(x, fill = b)) + geom_histogram(bins = 50, alpha = .5) +
  scale_fill_manual(values = c('red', 'blue')) + 
  guides(fill = guide_legend(override.aes = list(alpha = .5))) + 
  scale_alpha_manual(values = c(.5), guide = FALSE) +
  labs(x = 'Predicted value', y = 'Frequency', title = 'Cross-validated predictions')
cv_plot

dfpred <- data.frame(x = c(model_best$pred), b = as.factor(b))
pred_plot <- ggplot(dfpred, aes(x, fill = b)) + geom_histogram(bins = 50, alpha = .5) +
  scale_fill_manual(values = c('red', 'blue')) + 
  guides(fill = guide_legend(override.aes = list(alpha = .5))) + 
  scale_alpha_manual(values = c(.5), guide = FALSE) +
  labs(x = 'Predicted value', y = 'Frequency', title = 'Fully-fit predictions')
pred_plot

grid.arrange(pred_plot, cv_plot, ncol = 2)


stargazer(cbind(names(model_best$M) , round(model_best$beta, 3)))
model_best$beta

set.seed(14)
spikes2 <- spikes[sample(1:64, 10)]
system.time(results_2 <- opt_bh(b, spikes2, 0, 2, h_len_max = 7))
system.time(results_3 <- opt_bh(b, spikes2, .75, 1.25, h_len_max = 7))
system.time(results_4 <- opt_bh(b, spikes2, .95, 1.05, h_len_max = 7))
system.time(results_5 <- opt_bh(b, spikes2, .05, .125, h_len_max = 7))

size <- 24 # 6 + 6 + 6 + 6
errors_2 <- list(size)
h_s_2 <- errors_2
ind <- 1
for(i in seq_along(results_2$models)){
  errors_2[[ind]] <- results_2$models[[i]]
  h_s_2[ind] <- results_2$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results_3$models)){
  errors_2[[ind]] <- results_3$models[[i]]
  h_s_2[ind] <- results_3$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results_4$models)){
  errors_2[[ind]] <- results_4$models[[i]]
  h_s_2[ind] <- results_4$h_s[[i]]
  ind <- ind + 1
}
for(i in seq_along(results_5$models)){
  errors_2[[ind]] <- results_5$models[[i]]
  h_s_2[ind] <- results_5$h_s[[i]]
  ind <- ind + 1
}

grid_model_plot(errors_2, h_s_2)

model_best2 <- errors[[which.min(h_s_2)]]

(pred_rate2 <- sum((model_best2$pred > .5) == b)/length(model_best2$pred))

n <- length(b)
cv3 <- rep(0, n)
guess3 <- cv
for(i in 1:n){
  M_h3 <- model_best2$M[-i,]
  new_model3 <- unlist(mca_thresh_optim(b[-i], M_s = list(M_h2)), recursive = FALSE)
  guess3[i] <- sigmoid(sum(new_model3$beta*model_best2$M[i,]))
  cv3[i] <- b[i] == (guess3[i] > .5)
  
}

(cv3_rate <- sum(cv3)/length(cv3))
dfcv3 <- data.frame(x = guess3, b = b)
cv3_plot <- ggplot(dfcv3, aes(x, fill = b)) + geom_histogram(bins = 50, alpha = .5) +
  scale_fill_manual(values = c('red', 'blue')) + 
  guides(fill = guide_legend(override.aes = list(alpha = .5))) + 
  scale_alpha_manual(values = c(.5), guide = FALSE) +
  labs(x = 'Predicted value', y = 'Frequency', title = 'Cross-validated predictions')


roots_lone <- process_roots(t_over_h(mu_h, c(1), h_best))
lone_h_best <- int_mu_h(roots_lone$a, roots_lone$b, c(1), h_best)

minfl <- c('AF1', 'C4', 'CPZ', 'F5', 'FC6', 'FCZ', 'FT8', 'FZ', 'OZ', 'P7', 'P8', 'PO7', 'Y')
most_infl <- c(-0.015, -0.02, -0.017, 0.011, 0.01, 0.01, -0.02, 0.011, -0.013, 0.013, 0.014, 0.011, 0.011) 
OR1 <- round(exp(lone_h_best*most_infl), 4) # gives odds ratios

ch_infl <- colMeans(model_best$M)
OR2 <- round(exp(ch_infl[names(ch_infl) %in% minfl]*most_infl), 4)

stargazer(cbind(most_infl, OR1, OR2))
