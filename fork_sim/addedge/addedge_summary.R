# Investigate added edges
library(igraph)
library(tidyverse)
library(ggplot2)
library(sp)
load('../fork_data.Rdata')

dir_path <- 'visgp_fit_addedge'
pattern <- 'edge_'
files <- list.files(path = dir_path, pattern = pattern, full.names = T)
edge_diffs <- sapply(files, readRDS) # numbers of edges added 

p <- ggplot(data.frame(edge_diffs), aes(x=edge_diffs)) + 
    geom_histogram(bins=25, fill="blue", color="black") + 
    ggtitle("Histogram of edges added by chordal\n completion for n=200 in the fork-shaped domain") 

ggsave("addedges.png", plot=p, width=10, height=6, dpi=300)


# maximum edge_diff is at id=242
G <- readRDS(file.path(dir_path, 'G_242.rds'))
G.new <- readRDS(file.path(dir_path, 'G.new_242.rds'))

# Get the edge lists
el_new <- get.edgelist(G.new)
el <- get.edgelist(G)

# Convert edge lists to data frames for easier manipulation
df_new <- as.data.frame(el_new, stringsAsFactors = FALSE)
df <- as.data.frame(el, stringsAsFactors = FALSE)

# Renaming the columns
names(df_new) <- c("V1", "V2")
names(df) <- c("V1", "V2")

# Find edges in G.new that are not in G
edges_in_new_not_in_G <- anti_join(df_new, df)

# Convert to integer
edges_in_new_not_in_G$V1 <- as.integer(edges_in_new_not_in_G$V1)
edges_in_new_not_in_G$V2 <- as.integer(edges_in_new_not_in_G$V2)

# Get the training points
set.seed(242)
n = 250
stdev = .1

ind.train = sample.int(10000, .8*n, replace = F)
noise.train = rnorm(.8*n, sd = stdev)

ind.test = sample.int(10000, .2*n, replace = F)
noise.test = rnorm(.2*n, sd = stdev)

grid.train = grid.train[ind.train,]


png("addedges_largest.png", width=800, height=600)

plot(water)
points(grid.train)

for(i in 1:nrow(edges_in_new_not_in_G)) {
  lines(grid.train[as.numeric(edges_in_new_not_in_G[i,1:2]),], col='green')
}

title(main="Largest quantity of edges added by \nchordal completion")

dev.off()
