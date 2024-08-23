# Calculate log likelihood for data under visGP model
# param: Vector of model parameters (range, nugget, sigma, mu)
# nu: Smoothness parameter for Matérn covariance
# y: Vector of observed data
# dists: Matrix of distances between locations
# cliques: List of cliques in the visibility graph
# seps: List of separators in the visibility graph
# coords: Matrix of spatial coordinates
cov_select_loglik = function(param, nu = nu, y, dists, cliques,
                             seps, coords) {
  range = exp(param[1]) %>% as.numeric()
  nugget = exp(param[2]) %>% as.numeric()
  sigma = exp(param[3]) %>% as.numeric()
  mu = param[4]
  names(range) = 'range'
  names(nugget) = 'nugget'
  names(sigma) = 'sigma'
  names(mu) = 'mean'
  print(c(range, nugget, sigma, mu))
  
  out = -Reduce('+', lapply(cliques, function(clique) {
    if (length(clique) == 0) return(0)
    cov_mat = geostatsp::matern(dists[clique, clique], param = c(range=as.numeric(range), variance=as.numeric(sigma)^2, shape = nu)) + 
      diag(nugget, length(clique))
    return(mvtnorm::dmvnorm(y[clique], mean = rep(mu, length(clique)), sigma = cov_mat, log=T))
  })) + Reduce('+', lapply(seps, function(sep) {
    if (length(sep) == 0) return(0)
    cov_mat = geostatsp::matern(dists[sep, sep], param = c(range=as.numeric(range), variance=as.numeric(sigma)^2, shape = nu)) + 
      diag(nugget, length(sep))
    return(mvtnorm::dmvnorm(y[sep], mean = rep(mu, length(sep)), sigma = cov_mat, log=T))
  }))
  
  if (is.nan(out) | out == Inf | out == -Inf) return(9e99)
  return(out)
}

# Given fitted visGP model, make prediction at new points
# new_loc: Coordinates of new location for prediction
# grid: Matrix of observations
# y: Vector of observed data values
# param: Vector of model parameters
# nu: Smoothness parameter for Matérn covariance
# n_neighbors: Number of neighbors to use in prediction
# A: Adjacency matrix of the visibility graph
# A_vec: Optional vector of adjacency information for new location
# D: Matrix of distances between locations
# water: SpatialPolygons object representing the water body
# thresh: Distance threshold for considering points as neighbors
# method: Method for calculation/neighbor selection ('euclidean', 'nearest_clique', etc.)
pred_neighbor = function(new_loc, grid, y, param, nu = .5, n_neighbors, A, A_vec = NULL, D, water, thresh = Inf, method = NULL) {
  require(igraph)
  require(fields)
  
  range = exp(param[1]) %>% as.numeric()
  nugget = exp(param[2]) %>% as.numeric()
  sigma = exp(param[3]) %>% as.numeric()
  mu = param[4]
  
  euc_dists = as.vector(fields::rdist(matrix(new_loc, nrow=1), grid))
  
  if (is.null(A_vec)) {
    require(rgeos)
    within_thresh = which(euc_dists < thresh)
    
    l <- vector("list", length(within_thresh))
    
    for (j in seq_along(within_thresh)) {
      l[[j]] <- Lines(list(Line(as.matrix(rbind(new_loc, grid[within_thresh[j],])))), as.character(within_thresh[j]))
    }
    
    con = rep(F, nrow(grid))
    con[within_thresh] = gCovers(water, SpatialLines(l), byid=T)
  } else {
    con = A_vec
  }
  
  adj_dists = euc_dists/con
  
  cutoff = sort(adj_dists)[n_neighbors]
  if (cutoff == Inf) {
    neighbor_ind = which(adj_dists < Inf)
  } else {
    neighbor_ind = which(adj_dists <= cutoff)[1:n_neighbors]
  }
  
  if (length(neighbor_ind) == 0) return(list(mu=NA, sd = NA))
  
  if (method == 'euclidean') {
    dist_nn = D[neighbor_ind, neighbor_ind]
    dist_nn = rbind(adj_dists[neighbor_ind], dist_nn)
    dist_nn = cbind(c(0, adj_dists[neighbor_ind]), dist_nn)
    
    L = nrow(dist_nn)
    
    cov_mat = geostatsp::matern(dist_nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov_mat[1,1]
    Sig12 = cov_mat[1, 2:L, drop=F]
    Sig22inv = solve(cov_mat[2:L, 2:L])
    
    mu_post = mu + Sig12 %*% Sig22inv %*% (y[neighbor_ind] - mu)
    sd_post = sqrt(Sig11 - Sig12 %*% Sig22inv %*% t(Sig12))
    
    return(list(mu = mu_post, sd = sd_post))
  }
  
  neighbor_adj = A[neighbor_ind, neighbor_ind]
  L = length(neighbor_ind)
  
  if (method == 'nearest_clique') {
    if (!all(neighbor_adj + diag(L) == 1)) {
      neighbor_ind_red = c()
      
      ord = neighbor_ind[sort(adj_dists[neighbor_ind], index.return=T)$ix]
      for (candidate in ord) {
        if (all(A[c(neighbor_ind_red, candidate), c(neighbor_ind_red, candidate)] +
                diag(1+length(neighbor_ind_red)) == 1)) {
          neighbor_ind_red = c(neighbor_ind_red, candidate)
        }
      }
    } else {
      neighbor_ind_red = neighbor_ind
    }
    
    dist_nn = D[neighbor_ind_red, neighbor_ind_red]
    dist_nn = rbind(adj_dists[neighbor_ind_red], dist_nn)
    dist_nn = cbind(c(0, adj_dists[neighbor_ind_red]), dist_nn)
    
    L = nrow(dist_nn)
    
    cov_mat = geostatsp::matern(dist_nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov_mat[1,1]
    Sig12 = cov_mat[1, 2:L, drop=F]
    Sig22inv = solve(cov_mat[2:L, 2:L])
    
    mu_post = mu + Sig12 %*% Sig22inv %*% (y[neighbor_ind_red] - mu)
    sd_post = sqrt(Sig11 - Sig12 %*% Sig22inv %*% t(Sig12))
    
    return(list(mu = mu_post, sd = sd_post))
  }
  
  g = graph_from_adjacency_matrix(neighbor_adj, 'undirected')
  g = set.vertex.attribute(g, "name", value=1:length(neighbor_ind))
  
  if (method == 'maxprec') {
    cliques_list = max_cliques(g)
  } else if (method == 'precweighted') {
    i = 1
    cliques_list = list()
    while (vcount(g) > 0) {
      clique = attr(largest_cliques(g)[[1]], 'name')
      cliques_list[[i]] = as.numeric(clique)
      g = delete.vertices(g, clique)
      i = i + 1
    }
  }  
  
  out = lapply(cliques_list, function(clique) {
    neighbor_ind_clique = neighbor_ind[clique]
    
    dist_nn = D[neighbor_ind_clique, neighbor_ind_clique]
    dist_nn = rbind(adj_dists[neighbor_ind_clique], dist_nn)
    dist_nn = cbind(c(0, adj_dists[neighbor_ind_clique]), dist_nn)
    
    L = nrow(dist_nn)
    
    cov_mat = geostatsp::matern(dist_nn, param = c(range=range, variance=sigma^2, shape = nu)) + 
      diag(nugget, L)
    
    Sig11 = cov_mat[1,1]
    Sig12 = cov_mat[1, 2:L, drop=F]
    Sig22inv = solve(cov_mat[2:L, 2:L])
    
    mu_post = mu + Sig12 %*% Sig22inv %*% (y[neighbor_ind_clique] - mu)
    prec_post = 1/(Sig11 - Sig12 %*% Sig22inv %*% t(Sig12))
    
    return(list(mu_post=mu_post, prec_post=prec_post))
  })
  
  mu_posts = sapply(out, function(x) x$mu_post)
  prec_posts = sapply(out, function(x) x$prec_post)
  
  if (method == 'precweighted') return(list(mu=sum(mu_posts*prec_posts)/sum(prec_posts), sd = (sum(prec_posts))^(-.5)))
  else if (method == 'maxprec') return(list(mu=mu_posts[which.max(prec_posts)], sd = (max(prec_posts))^(-.5)))
}

# Fit visGP model by maximum likelihood
# coords: Matrix of spatial coordinates
# adj: Adjacency matrix of the visibility graph
# y: Vector of observed data
# nu: Smoothness parameter for Matérn covariance
fit_water = function(coords, adj, y, nu = .5) {
  require(igraph)
  require(gRbase)
  require(BRISC)
  
  D = as.matrix(dist(coords))
  
  g_unchord = graph_from_adjacency_matrix(adj, 'undirected')
  V(g_unchord)$name = 1:(length(y))
  g_chord = is_chordal(g_unchord, newgraph=T)$newgraph
  
  rip = rip(g_chord)
  
  # avoid numerical issues if simulation run has dense locations
  if (exists('dense') && isTRUE(dense) && length(y > 2500)) {
    range_start <- 1
    nugget_start <- 1
    sigma_start <- 1
    mu_start <- 0
  } else {
    tryCatch({
      BRISC_fit = BRISC_estimation(y=y, coords=coords)
      range_start = 1/BRISC_fit$Theta[3]
      nugget_start = BRISC_fit$Theta[2] + .01
      sigma_start = sqrt(BRISC_fit$Theta[1])
      mu_start = BRISC_fit$Beta
    }, error = function(e) {
      warning("Error in BRISC fit, initializing parameters to default values.")
      range_start <<- 1
      nugget_start <<- 1
      sigma_start <<- 1
      mu_start <<- 0
    })
  }
  
  cliques = sapply(rip$cliques, as.numeric)
  seps = sapply(rip$separators, as.numeric)
  
  if (length(y) > 2500) {
    fit = sgd(est_start= c(log(c(range_start, nugget_start, sigma_start)), mu_start),
              y_clique = lapply(cliques, function(c) y[c]), y_sep = lapply(seps, function(s) y[s]), 
              D_clique = lapply(cliques, function(c) D[c,c]), D_sep = lapply(seps, function(s) D[s,s]), 
              alpha = .025, beta = .99, epsilon = 1e-9)
  } else {
    fit = optim(par = c(log(c(range_start,nugget_start,sigma_start)), mu_start), fn = cov_select_loglik, nu = nu, y = y, dists = as.matrix(stats::dist(coords)), cliques=cliques,
                seps=seps, coords = coords)
  }
  
  return(fit)
}

# Calculate covariance matrix for given parameters
# covpar: Vector of covariance parameters
# D: Matrix of distances between locations
Sigma_fun = function(covpar, D) {
  exp(2*covpar[3])*exp(-D*2*exp(-covpar[1])) + 
    exp(covpar[2])*diag(nrow(D))
}

# Calculate derivative of covariance matrix w.r.t. parameters
# covpar: Vector of covariance parameters
# D: Matrix of distances between locations 
dSigma_fun = function(covpar, D) {
  c1 = 2 * exp(-covpar[1])
  c2 = exp(2 * covpar[3])
  c3 = exp(-D * c1)
  c4 = exp(covpar[2])
  list(c2 * c3 * D * c1, c4 * diag(nrow(D)), c2 * c3 * 2 * c2)
}

# Calculate gradient of log-likelihood
# y: Vector of observed data
# covpar: Vector of covariance parameters
# mu: Mean parameter
# D: Matrix of distances between locations
grad = function(y, covpar, mu, D) {
  if (length(y) == 0) return(rep(0, length(covpar)+1))
  Sigma = Sigma_fun(covpar, D)
  W = Matrix::chol2inv(Matrix::chol(Sigma))
  dS = dSigma_fun(covpar, D)
  dW = lapply(dS, function(x) -W %*% x %*% W)
  ydiff = y - mu
  grad_cov = sapply(1:3, function(i) -1/2*( t(ydiff) %*% dW[[i]] %*% (ydiff)) -
                      1/2*sum( diag(W %*% dS[[i]]) ) )
  grad_mu = sum(W %*% (ydiff))
  return(c(grad_cov, grad_mu))
}

# Perform stochastic gradient descent optimization
# est_start: Initial parameter estimates
# y_clique, y_sep: Observed data for cliques and separators
# D_clique, D_sep: Distance matrices for cliques and separators
# alpha: Learning rate
# beta: Momentum parameter
# epsilon: Small constant for numerical stability
# count_max: Maximum number of iterations
sgd = function(est_start, y_clique, y_sep, D_clique, D_sep, alpha = .025, beta = .99, epsilon = 1e-9, count_max = 5000) {
  m = length(y_clique)
  est_cur = est_start
  v_t <- numeric(4)
  count = 0
  repeat {
    ord <- sample.int(m)
    
    for (i in ord) {
      count = count + 1
      g <- grad(y_clique[[i]], est_cur[1:3], est_cur[4], D_clique[[i]]) -
        grad(y_sep[[i]], est_cur[1:3], est_cur[4], D_sep[[i]])
      
      v_t <- beta * v_t + (1 - beta) * (g * g)       # Update second moment estimate
      step <- alpha * g / (sqrt(v_t) + epsilon)      # Update step
      est_new <- est_cur + step
      
      print(c(exp(c(est_new)[1:3]), est_new[4]))
      est_cur <- est_new
      if (count >= count_max) break
    }
    return(list(par=est_cur))
  }
}