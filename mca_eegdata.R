library(eegkitdata)
library(ggplot2)
# data(eegdata)
setwd('/Users/adamdavis/Documents/Senior/Thesis/Code/')
source('mca_functions.R')

train_wd <- '/Users/adamdavis/Documents/Senior/Thesis/Code/EEGData/SMNI_CMI_TRAIN/'
test_wd <- '/Users/adamdavis/Documents/Senior/Thesis/Code/EEGData/SMNI_CMI_TEST/'
full_wd <- '/Users/adamdavis/Documents/Senior/Thesis/Code/EEGData/eeg_full/'

train_eeg <- geteegdata(train_wd, cond = 'S1')
test_eeg <- geteegdata(test_wd, cond = 'S1')
eeg_full <- geteegdata(full_wd, cond = 'S1')

eegdata <- within(eegdata, {
  unique_id <- paste(subject, trial, sep = '_')
  b <- group == 'a'
})

first_try <- eegdata[eegdata$subject == unique(eegdata$subject)[1] & 
                       eegdata$trial == 0 &
                       eegdata$channel == 'FP1',]
first_try


spik <- id_spikes(first_try)

spik

plot(first_try$time, first_try$voltage)
lines(first_try$time, first_try$voltage)
points(unname(unlist(spik))-1, rep(min(first_try$voltage[unname(unlist(spik))]), length(unname(unlist(spik)))), col = 'red')



first_try$voltage[unlist(spik)]

unique(eegdata$subject)

all_spikes <- id_spikes(eegdata)
all_2 <- all_spikes %>% transpose()
b <- aggregate(b ~ subject + trial, data = eegdata, FUN = mean)
b <- b[,3]

h_s <- seq(0, 1, .1)
M_s <- lapply(h_s, function(x){
  return(make_M(all_2, x))
})

model_1 <- mca_thresh_optim(b, h_s = h_s, M_s = M_s)
grid_model_plot(model_1, h_s)

# without the final thing bc interprettion is fuzzy:
ind_for_f <- seq(3, 183, 3)
M_s2 <- lapply(M_s, function(x){
  return(x[,-ind_for_f])
})

M_s2

model_2 <- mca_thresh_optim(b, h_s = h_s, M_s = M_s2)
grid_model_plot(model_2, h_s)

sum((model_2[[1]]$pred > .5) == b)

new_M <- make_M(all_2, .05)
f