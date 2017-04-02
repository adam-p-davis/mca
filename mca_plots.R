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
mu_h2 <- function(t){
  return(mu_h(t, t_starts = c(0, 2, 3), h = 0))
}

n <- 10000
super_plot <- data.frame(x = seq(-1, 8, length.out = n))
super_plot <- within(super_plot, {
  mu_h <- mu_h2(x)
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
