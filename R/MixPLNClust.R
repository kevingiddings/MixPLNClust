#' Clustering method for discrete count data with missing entries.
#'
#' This function takes in a matrix of partially observed count data, as well as a number of clusters \eqn{G}, and attempts to cluster
#' the observed data into \eqn{G} clusters according to a mixture of multivariate Poisson log-normal distributions.
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
#' \item ARI - Adjusted Rand index of estimated clustering against true cluster assignment, if true cluster assignment is provided. 
#' \item z - Matrix of percent chance for each observation to be in the respective cluster.
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
#'    PLNClust_results[[cluster]] <- MixPLNClust(MixPLN_data$Counts,cluster)
#' }
#' mclust::map(PLNClust_results[[2]]$z)
#' @export
MixPLNClust <- function(count_matrix,G,group_labels=NULL,parameter_space=NULL,
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

  while (checks==0){
    for (g in 1:G){
      GX[[g]] <- list()
      dGX[[g]] <- list()
      z_S[[g]] <- list()

      ## update approximation q parameters per cluster
      tryCatch({
        g_params <- update_g_params(GX[[g]], dGX[[g]], z_S[[g]], m[[g]], O_list,
        start[[g]], S[[g]], iOsigO[[g]], lib_mat, ini_Y, mu[[g]], z, g, N)
      }, error = function(e) {
        print(paste("Error:", e$message))
      })

      S[[g]] <- g_params$S_g
      dGX[[g]] <- g_params$dGX_g
      z_S[[g]] <- g_params$z_S_g
      GX[[g]] <- g_params$GX_g
      start[[g]] <- g_params$m_g

      # updating mu
      tryCatch({
        mu_new <- update_mu(m[[g]], O_list, iOsigO[[g]], mu[[g]], z, d, g, N)
      }, error = function(e) {
        print(paste("Error:", e$message))
      })
      mu[[g]] <- mu_new$mu_g


      # updating sigma
      tryCatch({
        sigma_new_g <- update_sig(m[[g]], sigma_new, sigma[[g]], O_list, iOsigO[[g]], mu[[g]], S[[g]], z, d, g, N, step_size)
      }, error = function(e) {
        print(paste("Error:", e$message))
      })
      sigma_new[[g]] <- sigma_new_g$sigma_new_g
      
    }

    ## SPD check. Update sigma_g,t to sigma_g,(t+1) if new sigma meets requirements
    tryCatch({
      PD_Check <- PD_check(sigma_new, G)
      if (PD_Check){
        sigma <- sigma_new
      }
      }, error = function(e) {
      print(paste("Error:", e$message))
    })

    ## calculate inverse sigma matrices
    invert_matrices(sigma, isigma, G)

    ## update percent of observations in each cluster
    pi_g <- compute_pi_g(z, N)
    
    ## update approximating q covariance param for each observation
    compute_iOsigO(O_list, sigma, iOsigO)
    
    
    tryCatch({
      results <- compute_F_matrices(S, m, O_list, mu, iOsigO, O_mat, ini_Y, pi_g)
      F_raw <- results$F_raw
      F <- results$F
      loglik[it] <- results$loglik
      z <- results$z
                }, error = function(e) {
      print(paste("Error:", e$message))
    })

    # aitkens acceleration check
    tryCatch({
      if (it > 3){
        aitkens <- aitkens_accel(it, loglik, aloglik)
        
        checks <- aitkens$checks
        exit_code <- aitkens$exit_code
        aloglik <- aitkens$a_loglik  
      }
    }, error = function(e) {
      print(paste("Error:", e$message))
    })
    
    
    #print(it)
    if (it>=max_iter){
      checks <- 1
      exit_code <- "Max iterations reached"
    }
    
    it <- it+1
  }
  k <- (G-1) +G*d +G*d*(d+1)/2
  BIC <- 2*loglik[length(loglik)] - k*log(N)
  ARI <- NULL
  if(!is.null(group_labels)){
    ARI <- mclust::adjustedRandIndex(mclust::map(z),group_labels)
  }
  ptm <- proc.time() - ptm
  Y[which(Y==-999)] <- NA
  return(list(Y=Y,pi_g=pi_g,mu=mu,sigma=sigma,z=z,ARI=ARI,loglik=loglik,BIC=BIC,kmeans=k_means,true=true,time=ptm,exit_code=exit_code))
}
