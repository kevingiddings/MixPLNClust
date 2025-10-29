#' Clustering method for discrete count data with missing entries.
#'
#' This function takes in a matrix of partially observed count data, as well as a number of clusters \eqn{G}, and attempts to cluster
#' the observed data into \eqn{G} clusters according to a mixture of multivariate Poisson log-normal distributions.
#' WARNING: This function performs all calculations in R, which is much slower than the version written with Rcpp (\code{\link[MixPLNClust:MixPLNClust]{MixPLNClust}}).
#' @param count_matrix A matrix of count data, potentially with missing values. Missing values should be denoted with NA values. Rows will be clustered. Rows with no observed values will be discarded.
#' @param G Positive integer specifying the number of clusters to test.
#' @param group_labels Optional. A vector of integers identifying the true cluster each row belongs to.
#' @param parameter_space Optional. Matrix of underlying parameter space that data was generated from.
#' @param gr_method Deprecated. Previously used to determine method of gradient descent. Defaults to using constant step size (CSS).
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
MixPLNClust_OLD <- function(count_matrix,G,group_labels=NULL,parameter_space=NULL,
                        gr_method="CSS",max_iter=1000, step_size = 0.0005,
                        calc_norm_factors = FALSE, custom_lib_mat = NULL,
                        init = "kmmeans"){
  
  if(is.null(count_matrix)){
    stop("Argument count_matrix is missing, with no default.")
  }else if(!is.matrix(count_matrix)){
    stop("Count data, count_matrix, must be a matrix.")
  }
  if(is.null(G)){
    stop("Argument G is missing, with no default.")
  }
  else if(G%%1 != 0 | G < 1){
    stop("Number of clusters, G, must be a positive integer.")
  }
  if(step_size <= 0 | !is.numeric(step_size)){
    stop("Step size, step_size, must be a positive number.")
  }
  if(max_iter <= 0 | !is.numeric(max_iter)){
    stop("Maximum number of iteration. max_iter, must be a positive number.")
  }
  ptm <- proc.time()
  ###### Parameter Updates ####

  Y <- count_matrix
  true <- group_labels
  true_par <- parameter_space

  exit_code <- "UNKNOWN STOPPING; POSSIBLE BUG"

  N <- nrow(Y) #sample size
  d <- ncol(Y) #dimension of the data

  #Create observation matrix
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
  


  #updated count matrix with entries determined to be missing set to 0.
  Y[is.na(O_mat)] <- NA

  ###Removes any observation with all values missing
  all_miss <- apply(Y, 1, function(x) all(is.na(x)))
  Y <- Y[ !all_miss, ]

  N <- nrow(Y) #sample size
  d <- ncol(Y) #dimension of the data

  O_list <- list()
  for (i in 1:N){
    O_list[[i]] <- diag(d)[which(O_mat[i,]==1),,drop=FALSE]
  }

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
  kept <- as.numeric(rownames(stats::na.omit(ini_Y)))
  k_means <- NULL
  if(init == "kmmeans"){
    k_means <- kmmeans::kmmeans(as.data.frame(log(ini_Y+1)),K=G,n.init = 100)[["partition"]]  ##Using k missing means to start to algorithm
    z <- mclust::unmap(k_means) ##Initial value for Z
  }
  if(init == "rand"){
    z_vec <- sample(c(1:G),N,replace=TRUE)
    z <- mclust::unmap(z_vec)
  }
  pi_g <- colSums(z)/N

  ###Initial value for Mu and Sigma
  ini_Y2 <- na.omit(ini_Y)
  ini_Y[is.na(ini_Y)] <- 0

  ###Initial value for Mu and Sigma
  for (g in 1:G){
    obs <- which(z[kept,g]==1)
    mu[[g]] <- colMeans(log(ini_Y2[obs,]+1/6))
    sigma[[g]] <- stats::var(log(ini_Y2[obs,]+1/6))
    isigma[[g]] <- MASS::ginv(sigma[[g]],tol=1e-20)
  }




  ###Initial value for m and S
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
      
      GX[[g]]<-list()
      dGX[[g]]<-list()
      z_S[[g]]<-list()
      for (i in 1:N){
        do<-nrow(O_list[[i]])
        dGX[[g]][[i]]<-diag(c(exp(O_list[[i]]%*%log(lib_mat)+start[[g]][[i]])+0.5*diag(S[[g]][[i]])),do)+iOsigO[[g]][[i]]
        S[[g]][[i]]<-solve(dGX[[g]][[i]]) ###S is updated here
        z_S[[g]][[i]]<-z[i,g]*S[[g]][[i]]
        GX[[g]][[i]]<-O_list[[i]]%*%ini_Y[i,]-exp(start[[g]][[i]]+O_list[[i]]%*%log(lib_mat)+0.5*diag(S[[g]][[i]]))-iOsigO[[g]][[i]]%*%(start[[g]][[i]]-O_list[[i]]%*%mu[[g]])
        m[[g]][[i]]<-start[[g]][[i]]+S[[g]][[i]]%*% GX[[g]][[i]] #m is updated here
      }
      start[[g]]<-m[[g]]
      
      
      for_mu_sig<-function(i){
        temp<-z[i,g]*t(O_list[[i]])%*%iOsigO[[g]][[i]]### This goes to the first part of sigma too
        num<-temp%*%m[[g]][[i]]
        den<-temp%*%O_list[[i]]
        ##binding such that first d by d are den and the last column is num
        return(cbind(den,num))
      }
      temp3<-Reduce("+",lapply(1:N,for_mu_sig))
      mu[[g]]<-c(solve(temp3[1:d,1:d])%*%temp3[,(d+1),drop=FALSE])
      
      
      ####Updating Sample covariance
      for_sig<-function(i){
        omega<-(m[[g]][[i]]-O_list[[i]]%*%mu[[g]])%*%t(m[[g]][[i]]-O_list[[i]]%*%mu[[g]])+S[[g]][[i]]
        forsig1<-z[i,g]*t(O_list[[i]])%*%iOsigO[[g]][[i]]%*%omega%*%iOsigO[[g]][[i]]%*%O_list[[i]]
        forsig2<-z[i,g]*t(O_list[[i]])%*%iOsigO[[g]][[i]]%*%O_list[[i]]
        return(-(forsig1-forsig2))
      }
      
      gr<-Reduce("+",lapply(1:N,for_sig))
      
      sigma_new[[g]]<-sigma[[g]]-step_size*gr
      
    }
    
    for_check<-unlist(lapply(sigma_new,function(x){eigen(x)$values}))
    
    if(all(for_check>0)){
      for(g in 1:G){
        sigma[[g]] <- sigma_new[[g]]
      }
    }
    
    for(g in 1:G){
      isigma[[g]] <- solve(sigma[[g]])
    }
    
    pi_g<-colSums(z)/N
    lib_mat_full<-matrix(rep(lib_mat,each=N),nrow=N) ###Matrix containing normaization factor so it makes easy to work with later.
    
    for (g in 1:G){
      for (i in 1:N){
        iOsigO[[g]][[i]]<-solve(as.matrix(O_list[[i]])%*%as.matrix(sigma[[g]])%*%t(as.matrix(O_list[[i]])))
      }
    }
    
    F<-matrix(NA,ncol=G,nrow=N)
    F_raw <- matrix(NA,ncol=G,nrow=N)
    for (g in 1:G){
      for (i in 1:N){
        F_raw[i,g] <- 0.5*log(det(S[[g]][[i]]))-0.5*t(m[[g]][[i]]-O_list[[i]]%*%mu[[g]])%*%iOsigO[[g]][[i]]%*%(m[[g]][[i]]-O_list[[i]]%*%mu[[g]])-sum(diag(iOsigO[[g]][[i]]%*%S[[g]][[i]]))+0.5*log(det(iOsigO[[g]][[i]]))+0.5*sum(na.omit(O_mat[i,]))+t(m[[g]][[i]])%*%O_list[[i]]%*%ini_Y[i,]-sum(exp(m[[g]][[i]]+0.5*diag(S[[g]][[i]]))+lfactorial(O_list[[i]]%*%ini_Y[i,]))
        F[i,g]<-pi_g[g]*exp(F_raw[i,g])
      }
    }
    # Find smallest positive value in F
    smallest_in_F <- sort(unique(F))[sort(unique(F))>0][1]
    #shift all values up by this smallest value to prevent division by zero
    F <- F + smallest_in_F
    
    loglik[it]<-sum(log(rowSums(F)))
    z<-F/rowSums(F)
    
    #### Numerical stability version. Seems to cause issue with logliklihood -> BIC -> selecting correct number of clusters
    #Fmax<-floor(apply(F_raw,1,max))
    #ll<<-sum(log(rowSums(exp(F_raw-Fmax)))+Fmax)
    #post<-exp(F_raw-Fmax-log(rowSums(exp(F_raw-Fmax))))
    #loglik[it]<-ll
    #z<-post
    
    if (it>3){
      #Aitkaine's stopping criterion
      if((loglik[it-1]-loglik[it-2])==0){
        checks<-1 
        exit_code <- "Log Likelihood equal for two iterations"
      }else{
        a<-(loglik[it]-loglik[it-1])/(loglik[it-1]-loglik[it-2])
        add_to<-(1/(1-a)*(loglik[it]-loglik[it-1]))
        aloglik[it]<-loglik[it-1]+add_to
        if(abs(aloglik[it]-loglik[it-1])<0.01){
          checks<-1
          exit_code <- "Aitken's acceleration converged"
        } 
      }
    }	
    print(it)
    it<-it+1
    if (it>=max_iter){
      checks<-1 
      exit_code <- "Max iterations reached"
    }   
  }
  k<- (G-1) +G*d +G*d*(d+1)/2
  BIC <- 2*loglik[length(loglik)] - k*log(N)   
  plot(loglik,type="l")
  ptm<-proc.time() - ptm
  Y[which(Y==-999)]<-NA
  return(list(Y=Y,pi_g=pi_g,mu=mu,sigma=sigma,z=z,loglik=loglik,BIC=BIC,kmeans=k_means,true=true,time=ptm,exit_code=exit_code))
}