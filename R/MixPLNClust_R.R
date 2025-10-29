#' Clustering method for discrete count data with missing entries.
#'
#' This function takes in a matrix of partially observed count data, as well as a number of clusters \eqn{G}, and attempts to cluster
#' the observed data into \eqn{G} clusters according to a mixture of multivariate Poisson log-normal distributions.
#' WARNING: This function performs all calculations in R, which is much slower than the version written with Rcpp (\code{\link[MixPLNClust:MixPLNClust]{MixPLNClust}}).
#' @param count_matrix A matrix of count data, potentially with missing values. Missing values should be denoted with NA values. Rows will be clustered. Rows with no observed values will be discarded.
#' @param G Positive integer specifying the number of clusters to test.
#' @param group_labels Optional. A vector of integers identifying the true cluster each row belongs to.
#' @param parameter_space Optional. Matrix of underlying parameter space that data was generated from.
#' @param max_iter Positive double. Max number of iterations to test. Default maximum for iterations is \eqn{1000}.
#' @param step_size Positive double. Step size for gradient descent method. Default step size is \eqn{0.0005}
#' @param calc_norm_factors Boolean. Should the data be scaled to account for different sampling depths?
#' @param custom_lib_mat Optional. Numerical vector of length d. Used in place of calc_norm_factors if user wants to specify a predefined scaling vector.
#' @param init Method for cluster initialization. Default is "kmmeans" to use \eqn{k}-missing-means. Other option is "rand", which randomly assigns all rows to 1 of the \eqn{G} groups with equal probability. Useful to prevent potential computational issues with large number of clusters. 
#' @return A list containing: 
#' \itemize{
#' \item Y - Count matrix supplied to MixPLNClust. 
#' \item pi_g - Vector of percent of observations in each cluster. 
#' \item mu - list of estimated means for each cluster. 
#' \item sigma - list of estimated covariance matrices for each cluster. 
#' \item z - Matrix of percent chance for each observation to be in the respective cluster.
#' \item ARI - Adjusted Rand index of estimated clustering against true cluster assignment, if true cluster assignment is provided.
#' \item loglik - Vector log likelihoods for each iteration. 
#' \item BIC - Bayesian inference criterion for the final iteration. 
#' \item kmeans - initial estimated clustering assignment according to \code{\link[kmmeans:kmmeans]{k missing means}}.
#' \item true - true cluster assignment (if supplied), 
#' \item time - time taken information 
#' \item exit_code - stopping reason (max iterations reached, log likelihood convergence, etc.)
#' }
#' @examples
#' #load in example data set
#' data(MixPLN_data)
#' 
#' clusters_to_test <- c(1,2,3)
#' PLNClust_results <- list()
#' 
#' #cluster with default settings
#' for(cluster in clusters_to_test){
#'    PLNClust_results[[cluster]] <- MixPLNClust_R(MixPLN_data$Counts,cluster)
#' }
#' mclust::map(PLNClust_results[[2]]$z)
#' @export
MixPLNClust_R <- function(count_matrix,G,group_labels=NULL,parameter_space=NULL,
                          max_iter=1000, step_size = 0.0005,
                          calc_norm_factors = FALSE, custom_lib_mat = NULL,
                          init = "kmmeans"){
  
  if(is.null(count_matrix)){
    stop("Argument count_matrix is missing, with no default.")
  }else if(!is.matrix(count_matrix)){
    stop("Count data, count_matrix, must be a matrix.")
  }
  if(is.null(G)){
    stop("Argument G is missing, with no default.")
  }else if(G%%1 != 0 || G < 1){
    stop("Number of clusters, G, must be a positive integer.")
  }
  if(step_size <= 0 || !is.numeric(step_size)){
    stop("Step size, step_size, must be a positive number.")
  }
  if(max_iter <= 0 || !is.numeric(max_iter)){
    stop("Maximum number of iterations, max_iter, must be a positive number.")
  }
  ptm <- proc.time()
  
  Y <- count_matrix
  true <- group_labels
  true_par <- parameter_space
  
  exit_code <- "UNKNOWN STOPPING; POSSIBLE BUG"
  
  N <- nrow(Y) #sample size
  d <- ncol(Y) #dimension of the data
  
  if(!is.null(group_labels) && length(group_labels) != N){
    stop("Group labels, group_labels, must be a vector of length N.")
  }
  
  #Create observation matrix
  #1 indicates an observation is present, 0 indicates missing
  O_mat <- matrix(NA,nrow = N,ncol = d)
  for(col in 1:d){
    for (row in 1:N){
      if(!is.na(Y[row,col])){
        O_mat[row,col] <- 1
      }
    }
  }
  
  #Drop rows where all entries were missing
  all_miss <- apply(O_mat, 1, function(x) all(is.na(x)))
  if(any(all_miss)){
    warning("Some rows contain only missing values! Dropping these rows.")
  }
  O_mat <- O_mat[ !all_miss, ]
  
  #extraction matrices for each individual
  O_list <- list()
  for (i in 1:N){
    O_list[[i]] <- diag(d)[which(O_mat[i,]==1),,drop=FALSE]
  }
  
  #preprocessing for non-simulated data that needs to be scaled
  Y_2 <- Y
  Y_2[is.na(Y_2)] <- 0
  all_miss <- apply(Y_2, 1, function(x) all(is.na(x)))
  Y_2 <- Y_2[!all_miss, ]
  lib_mat <- rep(1, d)
  if (calc_norm_factors) {
    lib_mat <- edgeR::calcNormFactors(Y_2)
  }
  if(!is.null(custom_lib_mat)){
    lib_mat <- custom_lib_mat
  }
  
  #### Initialization ###
  mu <- list()
  psi <- list()
  lambda <- list()
  sigma <- list()
  isigma <- list()
  sigma_new <- list()
  m <- list()
  S <- list()
  P <- list()
  Q <- list()
  
  ###Other intermediate items initialized
  start <- list()
  Sk <- array(0, c(d,d,G) )
  GX <- list()
  dGX <- list()
  iOsigO <- list()
  
  z_S <- list()
  z_SO <- list()
  z_DO <- list()
  
  ini_Y <- Y
  ini_Y[ini_Y==-999] <- NA
  rownames(ini_Y) <- 1:nrow(ini_Y)
  #kept = rows with complete data
  kept <- as.numeric(rownames(stats::na.omit(ini_Y)))
  k_means <- NULL
  if(init == "kmmeans"){
    k_means <- kmmeans::kmmeans(as.data.frame(log(ini_Y+1)),K=G,n.init = 100)[["partition"]]  ##Using k missing means for starting cluster partition
    z <- mclust::unmap(k_means) ##Initial value for Z
  }
  if(init == "rand"){
    z_vec <- sample(c(1:G),N,replace=TRUE) #if random initialization selected, just randomly assign clusters
    z <- mclust::unmap(z_vec)
  }
  pi_g <- colSums(z)/N
  
  
  ini_Y2 <- na.omit(ini_Y)
  ini_Y[is.na(ini_Y)] <- 0
  
  ###Initial value for Mu and Sigma
  for (g in 1:G){
    obs <- which(z[kept,g]==1)
    mu[[g]] <- colMeans(log(ini_Y2[obs,]+1/6))
    sigma[[g]] <- stats::var(log(ini_Y2[obs,]+1/6))
    isigma[[g]] <- MASS::ginv(sigma[[g]],tol=1e-20) #inverse sigma
  }
  
  
  
  
  ###Initial value for m and S. These are parameters for approximating density q
  for (g in 1:G){
    S[[g]] <- list()
    start[[g]] <- list()
    m[[g]] <- list()
    iOsigO[[g]] <- list()
    for (i in 1:N){
      do <- nrow(O_list[[i]])
      start[[g]][[i]] <- log(O_list[[i]]%*%ini_Y[i,]+1/6) ###Starting value for M
      m[[g]][[i]] <- start[[g]][[i]]
      S[[g]][[i]] <- diag(do)*0.000000001
      iOsigO[[g]][[i]] <- MASS::ginv(as.matrix(O_list[[i]])%*%as.matrix(sigma[[g]])%*%t(as.matrix(O_list[[i]])),tol=1e-20)
    }
  }
  
  
  checks <- 0
  it <- 1
  aloglik <- NULL
  loglik <- NULL ##Log likelihood is stored in this vector to check for convergence.
  aloglik[1:3] <- 0
  
  
  ## ---------- helpers ----------
  
  #trace of matrix
  .tr <- function(M) {
    sum(diag(M))
  }
  #ensure symmetry in face of numerical imprecision. Keeps quadratic forms equal
  .symmetrize <- function(M) {
    0.5 * (M + t(M))
  }
  
  
  update_g_params <- function(GX_g, dGX_g, z_S_g, m_g, O_list, start_g, S_g, iOsigO_g,
                              lib_mat, Y, mu_g, z, g, N) {
    if (length(S_g) < N){
      S_g <- rep(list(matrix(,0,0)), N)
    }
    if (length(dGX_g) < N){
      dGX_g <- rep(list(matrix(,0,0)), N)
    }
    if (length(z_S_g) < N){
      z_S_g <- rep(list(matrix(,0,0)), N)
    }
    if (length(GX_g) < N){
      GX_g <- rep(list(numeric()), N)
    }
    if (length(m_g) < N){
      m_g <- rep(list(numeric()), N)
    }
    
    ## clamp to [1e-16, +inf)
    lib_mat <- pmax(as.numeric(lib_mat), 1e-16)
    
    
    for (i in seq_len(N)) {
      O_i <- O_list[[i]]
      start_i <- as.numeric(start_g[[i]])
      S_i <- S_g[[i]]
      iOsigO_i <- iOsigO_g[[i]]
      
      ## dGX = diag(exp(O_i %*% log(lib_mat) + start) + 0.5*diag(S)) + iOsigO
      x_lin <- as.numeric(O_i %*% log(lib_mat)) + start_i
      diag_vec <- exp(x_lin) + 0.5 * diag(S_i)
      dGX_i <- diag(diag_vec, nrow = length(diag_vec)) + iOsigO_i
      S_i <- solve(dGX_i)
      
      z_S_i <- as.numeric(z[i, g]) * S_i
      
      ## GX = O_i * Y_i - exp(start + O_i*log(lib) + 0.5*diag(S)) - iOsigO*(start - O_i*mu)
      Y_i <- as.numeric(Y[i, ])                # row i
      term1 <- as.numeric(O_i %*% Y_i)           # vector
      term2 <- exp(start_i + as.numeric(O_i %*% log(lib_mat)) + 0.5 * diag(S_i))
      term3 <- iOsigO_i %*% (start_i - as.numeric(O_i %*% mu_g))
      GX_i <- term1 - term2 - as.numeric(term3)
      
      m_i <- start_i + S_i %*% GX_i
      m_i <- as.numeric(m_i)
      
      S_g[[i]] <- S_i
      dGX_g[[i]] <- dGX_i
      z_S_g[[i]] <- z_S_i
      GX_g[[i]] <- GX_i
      m_g[[i]] <- m_i
    }
    
    list(S_g = S_g, dGX_g = dGX_g, z_S_g = z_S_g, GX_g = GX_g, m_g = m_g)
  }
  
  
  .for_mu <- function(i, g, iOsigO_g, z, O_list, m_g_i) {
    O_i <- O_list[[i]]
    iOi <- iOsigO_g[[i]]
    temp <- as.numeric(z[i, g]) * t(O_i) %*% iOi
    num <- temp %*% m_g_i
    den <- temp %*% O_i
    cbind(den, num)  # d x (d+1)
  }
  
  
  update_mu <- function(m_g, O_list, iOsigO_g, mu_g, z, d, g, N) {
    mu_temp <- matrix(0, nrow = d, ncol = d + 1)
    for (i in seq_len(N)) {
      mu_temp <- mu_temp + .for_mu(i, g, iOsigO_g, z, O_list, m_g[[i]])
    }
    den <- mu_temp[, 1:d, drop = FALSE]
    num <- mu_temp[, d + 1, drop = FALSE]
    mu_g <- solve(den, num)                # d x 1
    mu_g <- as.numeric(mu_g)
    list(mu_g = mu_g)
  }
  
  
  .for_sig <- function(i, g, iOsigO_g, z, O_list, m_g_i, mu_g, S_g_i) {
    O_i <- O_list[[i]]
    iOi <- iOsigO_g[[i]]
    diff <- m_g_i - as.numeric(O_i %*% mu_g)
    omega <- tcrossprod(diff, diff) + S_g_i
    for1 <- as.numeric(z[i, g]) * t(O_i) %*% iOi %*% omega %*% iOi %*% O_i
    for2 <- as.numeric(z[i, g]) * t(O_i) %*% iOi %*% O_i
    -(for1 - for2)
  }
  
  
  update_sig <- function(m_g, sigma_new, sigma_g, O_list, iOsigO_g, mu_g, S_g, z, d, g, N, step) {
    gr <- matrix(0, d, d)
    for (i in seq_len(N)) {
      mi <- as.numeric(m_g[[i]])
      Si <- S_g[[i]]
      gr <- gr + .for_sig(i, g, iOsigO_g, z, O_list, mi, mu_g, Si)
    }
    sigma_new_g <- sigma_g - step * gr
    list(sigma_new_g = sigma_new_g)
  }
  
  
  PD_check <- function(sigma_new, G) {
    for (g in seq_len(G)) {
      Sg <- sigma_new[[g]]
      Sg <- .symmetrize(Sg)
      ev <- eigen(Sg, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev <= 0) || any(!is.finite(ev))) return(FALSE)
    }
    TRUE
  }
  
  
  invert_matrices <- function(sigma, i_sigma, G) {
    for (g in seq_len(G)) {
      i_sigma[[g]] <- chol2inv(chol(sigma[[g]]))
    }
    invisible(NULL)
  }
  
  
  compute_pi_g <- function(z, N) {
    colSums(z) / as.numeric(N)
  }
  
  compute_iOsigO <- function(O_list, sigma, iOsigO) {
    G <- length(sigma)
    N <- length(O_list)
    for (g in seq_len(G)) {
      Sg <- sigma[[g]]
      iOsigO_g <- vector("list", N)
      for (i in seq_len(N)) {
        O_i <- O_list[[i]]
        mid <- O_i %*% Sg %*% t(O_i)
        iOsigO_g[[i]] <- chol2inv(chol(mid))
      }
      iOsigO[[g]] <- iOsigO_g
    }
    invisible(NULL)
  }
  
  
  compute_F_matrices <- function(S, m, O_list, mu, iOsigO, O_mat, Y, pi_g) {
    G <- length(S)
    N <- length(O_list)
    
    F_raw <- matrix(NaN, N, G)
    F <- matrix(NaN, N, G)
    
    log_pi_g <- log(pi_g)
    
    for (g in seq_len(G)) {
      S_g <- S[[g]]
      m_g <- m[[g]]
      iO_g <- iOsigO[[g]]
      mu_g <- as.numeric(mu[[g]])
      
      for (i in seq_len(N)) {
        ## count of non-NA in O_mat row i
        O_sum <- sum(!is.na(O_mat[i, ]))
        
        Sgi <- S_g[[i]]
        iOi <- iO_g[[i]]
        O_i <- O_list[[i]]
        mgi <- as.numeric(m_g[[i]])
        Yi <- as.numeric(Y[i, ])
        
        term_det <- 0.5 * (log(det(Sgi)) + log(det(iOi)))
        diff <- mgi - as.numeric(O_i %*% mu_g)
        quad <- 0.5 * as.numeric(t(diff) %*% iOi %*% diff)
        tr_term <- .tr(iOi %*% Sgi)
        part1 <- term_det - quad - tr_term + 0.5 * O_sum
        part2 <- as.numeric(t(mgi) %*% O_i %*% Yi)
        
        ## vector of counts for lgamma: O_i %*% Y_i  (k x 1)
        oy <- as.numeric(O_i %*% Yi)
        lg_term <- sum(lgamma(oy + 1))
        
        ## exp(m + .5*diag(S)) vector sum
        exp_term <- sum(exp(mgi + 0.5 * diag(Sgi)))
        
        F_raw[i, g] <- part1 + part2 - (exp_term + lg_term)
      }
    }
    
    ## softmax row-wise with stabilization
    for (i in seq_len(N)) {
      lp   <- F_raw[i, ] + log_pi_g
      mLP  <- max(lp)
      w    <- exp(lp - mLP)
      F[i, ] <- w / sum(w)
    }
    
    ## log-likelihood via rowwise log-sum-exp
    loglik <- 0
    for (i in seq_len(N)) {
      lp   <- F_raw[i, ] + log_pi_g
      mLP  <- max(lp)
      loglik <- loglik + (mLP + log(sum(exp(lp - mLP))))
    }
    
    z <- F
    list(F_raw = F_raw, F = F, loglik = loglik, z = z)
  }
  
  
  aitkens_accel <- function(it, loglik, a_loglik) {
    checks <- 0L
    exit_code <- ""
    it <- it - 1L
    
    if (it >= length(loglik))   length(loglik)  <- it + 1L
    if (it >= length(a_loglik)) length(a_loglik) <- it + 1L
    
    if (it > 3L) {
      if ((loglik[it - 1L] - loglik[it - 2L]) == 0) {
        checks <- 1L
        exit_code <- "Log Likelihood equal for two iterations"
      } else {
        a <- (loglik[it] - loglik[it - 1L]) / (loglik[it - 1L] - loglik[it - 2L])
        add_to <- (1L / (1L - a)) * (loglik[it] - loglik[it - 1L])
        a_loglik[it] <- loglik[it - 1L] + add_to
        if (abs(a_loglik[it] - loglik[it - 1L]) < 0.01) {
          checks <- 1L
          exit_code <- "Aitken's acceleration converged"
        }
      }
      return(list(checks = checks, exit_code = exit_code, a_loglik = a_loglik))
    } else {
      return(list(checks = 0L, exit_code = "Iteration index too small", a_loglik = a_loglik))
    }
  }
  
  
  
  while (checks == 0) {
    for (g in seq_len(G)) {
      GX[[g]] <- list()
      dGX[[g]] <- list()
      z_S[[g]] <- list()
      
      ## update approximation q parameters per cluster
      g_params <- tryCatch(
        update_g_params(GX[[g]], dGX[[g]], z_S[[g]], m[[g]], O_list,
                        start[[g]], S[[g]], iOsigO[[g]], lib_mat, ini_Y,
                        mu[[g]], z, g, N),
        error = function(e) { message("Error: ", e$message); NULL }
      )
      if (is.null(g_params)) break
      
      S[[g]] <- g_params$S_g
      dGX[[g]] <- g_params$dGX_g
      z_S[[g]] <- g_params$z_S_g
      GX[[g]] <- g_params$GX_g
      start[[g]] <- g_params$m_g
      
      ## update mu
      mu_new <- tryCatch(
        update_mu(m[[g]], O_list, iOsigO[[g]], mu[[g]], z, d, g, N),
        error = function(e) { message("Error: ", e$message); NULL }
      )
      if (is.null(mu_new)) break
      mu[[g]] <- mu_new$mu_g
      
      ## update sigma
      sigma_new_g <- tryCatch(
        update_sig(m[[g]], sigma_new, sigma[[g]], O_list, iOsigO[[g]], mu[[g]],
                   S[[g]], z, d, g, N, step_size),
        error = function(e) { message("Error: ", e$message); NULL }
      )
      if (is.null(sigma_new_g)) break
      sigma_new[[g]] <- sigma_new_g$sigma_new_g
    }
    
    ## SPD check. Update sigma_g,t to sigma_g,(t+1) if new sigma meets requirements
    PD_Check <- tryCatch(PD_check(sigma_new, G),
                         error = function(e) { message("Error: ", e$message); FALSE })
    if (isTRUE(PD_Check)) sigma <- sigma_new
    
    ## calculate inverse sigma matrices
    invert_matrices(sigma, isigma, G)
    
    ## update percent of observations in each cluster
    pi_g <- compute_pi_g(z, N)
    
    ## update approximating q covariance param for each observation
    compute_iOsigO(O_list, sigma, iOsigO)
    
    ## E-step like + loglik
    results <- tryCatch(
      compute_F_matrices(S, m, O_list, mu, iOsigO, O_mat, ini_Y, pi_g),
      error = function(e) { message("Error: ", e$message); NULL }
    )
    if (is.null(results)) break
    F_raw <- results$F_raw
    F <- results$F
    loglik[it] <- results$loglik
    z <- results$z
    
    ## Aitken’s acceleration
    if (it > 3) {
      aa <- tryCatch(
        aitkens_accel(it, loglik, aloglik),
        error = function(e) { message("Error: ", e$message); NULL }
      )
      if (!is.null(aa)) {
        checks <- aa$checks
        exit_code <- aa$exit_code
        aloglik <- aa$a_loglik
      }
    }
    
    if (it >= max_iter) {
      checks <- 1
      exit_code <- "Max iterations reached"
    }
    
    it <- it + 1
  }
  k <- (G-1) +G*d +G*d*(d+1)/2
  BIC <- 2*loglik[length(loglik)] - k*log(N)   
  
  ARI <- NULL
  if(!is.null(group_labels)){
    ARI <- mclust::adjustedRandIndex(mclust::map(z),group_labels)
  }
  
  ptm<-proc.time() - ptm
  Y[which(Y==-999)]<-NA
  return(list(Y=Y,pi_g=pi_g,mu=mu,sigma=sigma,z=z,ARI=ARI,loglik=loglik,BIC=BIC,kmeans=k_means,true=true,time=ptm,exit_code=exit_code))
}
