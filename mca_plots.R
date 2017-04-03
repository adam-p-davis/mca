library(ggplot2)
library(gridExtra)
library(plotmath)
setwd('/Users/adamdavis/Documents/Senior/Thesis/Code/')
source('mca_functions.R')


# single impulse response
ggplot(data.frame(x = c(-1, 5)), aes(x)) + 
  stat_function(fun = f, geom = 'line', n = 10000) + 
  labs(title = 'Single impulse at t = 0', 
       x = 't', y = 'Impulse response')


# integrated impulse response
mu_h3 <- function(t){
  return(mu_h(t, t_starts = c(0, 2, 3), h = 0))
}

n <- 10000
super_plot <- data.frame(x = seq(-1, 8, length.out = n))
super_plot <- within(super_plot, {
  mu_h <- mu_h3(x)
  h <- rep(.5, n)
  h2 <- rep(.75, n)
  h3 <- rep(.25, n)
  h4 <- rep(1, n)
})

times <- c(0, 2, 3)
roots1 <- process_roots(t_over_h(mu_h, times, .5))
v1 <- int_mu_h(roots1$a, roots1$b, times, .5)

roots2 <- process_roots(t_over_h(mu_h, times, .75))
v2 <- int_mu_h(roots2$a, roots2$b, times, .75)

roots3 <- process_roots(t_over_h(mu_h, times, .25))
v3 <- int_mu_h(roots3$a, roots3$b, times, .25)

roots4 <- process_roots(t_over_h(mu_h, times, 1))
v4 <- int_mu_h(roots4$a, roots4$b, times, 1)

h25 <- ggplot(super_plot, aes(x = x, y = mu_h)) + 
  labs(title = 'Superimposed impulses at t = 0, 2, 3; h = .25', 
       x = 't', y = 'Impulse response') +
  geom_line(aes(y = mu_h)) + geom_line(aes(y = h3), size = 1.2) + 
  geom_ribbon(data = subset(super_plot, mu_h > h3), aes(ymin = h3, ymax = mu_h), fill = 'blue', alpha = 0.5) +
  scale_fill_manual(values=c(clear,blue)) + 
  annotate('text', x = 5, y = .65, label = paste('M(h) =', round(v3,2)))

h5 <- ggplot(super_plot, aes(x = x, y = mu_h)) + 
  labs(title = 'h = .5', 
       x = 't', y = 'Impulse response') +
  geom_line(aes(y = mu_h)) + geom_line(aes(y = h), size = 1.2) + 
  geom_ribbon(data = subset(super_plot, mu_h > h), aes(ymin = h, ymax = mu_h), fill = 'blue', alpha = 0.5) +
  scale_fill_manual(values=c(clear,blue)) +
  annotate('text', x = 5, y = .65, label = paste('M(h) =', round(v1,2)))

h75 <- ggplot(super_plot, aes(x = x, y = mu_h)) + 
  labs(title = 'h = .75', 
       x = 't', y = 'Impulse response') +
  geom_line(aes(y = mu_h)) + geom_line(aes(y = h2), size = 1.2) + 
  geom_ribbon(data = subset(super_plot, mu_h > h2), aes(ymin = h2, ymax = mu_h), fill = 'blue', alpha = 0.5) +
  scale_fill_manual(values=c(clear,blue)) + 
  annotate('text', x = 5, y = .65, label = paste('M(h) =', round(v2,2)))

h1 <- ggplot(super_plot, aes(x = x, y = mu_h)) + 
  labs(title = 'h = 1', 
       x = 't', y = 'Impulse response') +
  geom_line(aes(y = mu_h)) + geom_line(aes(y = h4), size = 1.2) + 
  geom_ribbon(data = subset(super_plot, mu_h > h4), aes(ymin = h4, ymax = mu_h), fill = 'blue', alpha = 0.5) +
  scale_fill_manual(values=c(clear,blue)) + 
  annotate('text', x = 5, y = .65, label = paste('M(h) =', round(v4,2)))

grid.arrange(h25, h5, h75, h1, ncol = 2)

# sample set of observations
m <- 3
len <- 256
t_lists <- list(m)
plot_list <- list(m)
h <- .5
for(i in 1:m){
  t_lists[[i]] <- which(sample(0:5, len, replace = TRUE) == 1)
  roots <- process_roots(t_over_h(mu_h, t_lists[[i]], h))
  v <- int_mu_h(roots$a, roots$b, t_lists[[i]], h)
  coords <- data.frame(x = seq(0, len, length.out = n))
  coords <- within(coords, {
    mu_h <- mu_h(x, t_starts = t_lists[[i]], h = 0)
    h5 <- rep(h, n)
  })
  plot_list[[i]] <- ggplot(coords, aes(x = x, y = mu_h)) + 
    geom_line(aes(y = mu_h)) + geom_line(aes(y = h5), size = 1.2) +
    labs(title = paste(m, 'sample observations'),
         x = 't', y = 'Impulse response') + 
    geom_ribbon(data = subset(coords, mu_h > h5), aes(x = x, ymin = h5, ymax = mu_h), fill = 'blue', alpha = 0.5) +
    scale_fill_manual(values=c(clear,blue)) + 
    annotate('text', x = 257, y = .65, label = paste0('M', i, '(h) = ', round(v,2))) +
    annotate('text', x = 265, y = .65, label = paste0('b', i, ' = ', sample(0:1)))
}

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], ncol = 1)
plot