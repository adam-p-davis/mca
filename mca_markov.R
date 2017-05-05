library(markovchain)

sample_states <- c('START', 'FB', 'GDN', 'CONV', 'NULL')
sample_matrix <- matrix(data = c(0, 2/3, 1/3, 0, 0,
                                0, 1/5, 3/5, 0, 1/5,
                                0, 2/4, 0, 1/2, 0,
                                0, 0, 0, 1, 0,
                                0, 0, 0, 0, 1), byrow = TRUE, nrow = 5)
sample <- new('markovchain', states = sample_states, 
              transitionMatrix = sample_matrix,
              name = 'Sample')
print(sample)
show(sample)

library(diagram)

plot(sample, package = 'diagram', box.size = 0.025)

plot(sample)
