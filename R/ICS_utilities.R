#' ICS based on QR decomposition in case of full rank matrix
#'
#' This algorithm performs a stable ICS for the following pair of scatter matrices:
#' cov-covW and full rank data.
#' cov is the usual variance-covariance matrix and covW is a one-step M-estimator.
#'
#' @param Xt numeric nxp data matrix or dataframe.
#' @param alpha parameter of the one-step M-estimator. By default equals to 1.
#' @param cf consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
#'
ICS_QR_full_rank <- function(Xt, alpha = 1, cf = 1){
  # Details
  n <- nrow(Xt)
  p <- ncol(Xt)

  # Normalization by the mean -----
  Xt.c = sweep(Xt, 2, colMeans(Xt), "-")

  # Permutation by rows ---------------------------------------------------
  # decreasing by infinity norm:  absolute maximum
  norm_inf <- apply(Xt.c, 1, function(x) max(abs(x)))
  order_rows <- order(norm_inf, decreasing = TRUE)
  Xt_row_per <- Xt.c[order_rows,]

  # QR decomposition of Xt with column pivoting from LAPACK --------------------
  qr_Xt <- qr( 1/sqrt(n-1)*Xt_row_per, LAPACK = TRUE)

  # R should be pxp
  R <- qr.R(qr_Xt)

  # Q should be nxp
  Q <- qr.Q(qr_Xt)
  # Checks
  # D <- crossprod(Q)
  # diag(D) # should be only 1

  # computation of d_i
  d_i <- (n-1)*apply(Q, 1, function(x) sum(x^2))

  # Spectral decomposition (SEP) of S2y --------------------------------------------------------------

  # Computation of S2
  # with general one-step M-estimators
  S2_y <- cf*1/n*(n-1)*t(Q) %*% diag(d_i^alpha) %*% Q


  # SEP of S2_y
  SEP_S2_y <- eigen(S2_y, symmetric = TRUE)
  U2 <- SEP_S2_y$vectors

  # Eigenvectors by rows ----
  Rinv <- qr.solve(R)
  B <- t(U2) %*% t(Rinv)

  # Checks
  # S1_x <- cov(Xt_row_per)
  # norm_S1 <- B%*%S1_x[qr_Xt$pivot,qr_Xt$pivot]%*%t(B)
  # diag(norm_S1)

  # Scores -----
  # Be careful to reorder the rows and cols
  scores_Z <-  Xt_row_per[order(order_rows),qr_Xt$pivot]%*%t(B)


  return(list(eigenvalues = SEP_S2_y$values,
              B = B,
              Z = scores_Z,
              di = d_i^alpha,
              alpha = alpha))
}

#' ICS based on QR decomposition in case of not full rank matrix
#'
#' This algorithm performs a stable ICS for the following pair of scatter matrices:
#' cov-covW in the case the data are not full rak.
#' cov is the usual variance-covariance matrix and covW is a one-step M-estimator.
#' In this case, first a dimension reduction is performed before calling the `ICS_QR_full_rank ` function.
#' The determination of the rank $q$ of the data is made based on the QR decomposition of the matrix.
#' Two methods are implemented: if "first", the singular values of R are compared to the first one and
#' only the ones higher than $max(n,p) *tol_eps$ are kept. If "previous" the singular values are compared
#' to the previous one and we stop as soon as the criteria is not met.
#'
#' @param Xt numeric nxp data matrix or dataframe.
#' @param ref The reference to estimate the rank q: either "previous" or "first".
#' @param tol_eps Value of tolerance 10^(-8)
#' @param alpha parameter of the one-step M-estimator. By default equals to 1.
#' @param cf consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
#'
ICS_QR_not_full_rank <- function(Xt, ref = "previous",  tol_eps = 10^(-8),
                                 alpha = 1, cf = 1){
  # Details
  n <- nrow(Xt)
  p <- ncol(Xt)

  # Normalization by the mean -----
  Xt.c = sweep(Xt, 2, colMeans(Xt), "-")

  # Permutation by rows ---------------------------------------------------
  # decreasing by infinity norm:  absolute maximum
  norm_inf <- apply(Xt.c, 1, function(x) max(abs(x)))
  order_rows <- order(norm_inf, decreasing = TRUE)
  Xt_row_per <- Xt.c[order_rows,]

  # QR decomposition of Xt with column pivoting from LAPACK --------------------
  qr_Xt <- qr(1/sqrt(n-1)*Xt_row_per, LAPACK = TRUE)


  # Estimation of rank q -------------------------------------------------
  # R is nxp, but with only zero for rows > p
  # the diag of R is already in decreasing order and is a good approximation
  # of the rank of X.c
  # R should be pxp
  R <- qr.R(qr_Xt)
  r_all <- abs(diag(R))

  if (ref =="first"){
    r_ratios <-  r_all/r_all[1]
    q <- sum(r_ratios > max(n,p) *tol_eps)
  }else{
    r_ratios <- r_all[2:length(r_all)]/r_all[1:(length(r_all)-1)]
    q <- which(r_ratios < max(n,p) *tol_eps)[1]
  }
  q <- ifelse(is.na(q), length(r_all), q)
  

  # Q should be nxp but we are only interested in nxq
  Q1 <- qr.Q(qr_Xt)[,1:q]
  # D <- crossprod(Q1)
  # diag(D) # should be only 1


  # QR decomposition of Rt ---------------------------------------------------
  R_q <- R[1:q, ]
  qr_R <- qr(t(R_q), LAPACK = TRUE)
  Tau <- qr.Q(qr_R)[1:q, ]
  Omega1 <- qr.R(qr_R)[1:q, 1:q]

  # New X tilde -------------------------------------------------------------
  # permutation matrices
  # permutation of rows
  Pi2 = data.frame(model.matrix(~ . -1, data = data.frame(row=as.character(order_rows))))
  Pi2 <- Pi2[,order(as.numeric(substr(colnames(Pi2), start = 4, stop = nchar(colnames(Pi2)))))]
  colnames(Pi2) <- rownames(Xt)
  
  # permutation of cols
  Pi3 = data.frame(model.matrix(~ . -1, data = data.frame(col=as.character( qr_R$pivot))))
  Pi3 <- t(Pi3[,order(as.numeric(substr(colnames(Pi3), start = 4, stop = nchar(colnames(Pi3)))))])

  X_tilde <- sqrt(n-1)* Tau %*% t(Pi3) %*% t(Q1)

  Xt_tilde <- t(Pi2) %*% t(X_tilde)


  # ICS based on QR ---------------------------------------------------------

  res_ICS_QR <- ICS_QR_full_rank(Xt = Xt_tilde, alpha = alpha, cf = cf)

  return(list(
    eigenvalues = res_ICS_QR$eigenvalues,
    B = res_ICS_QR$B,
    Z = res_ICS_QR$Z,
    Xt_tilde =  Xt_tilde,
    q = q,
    qi = res_ICS_QR$qi,
    alpha = res_ICS_QR$alpha,
    r_all =  r_all))
}




#'  One-step M-estimator covW
#'
#' @param X numeric nxp data matrix or dataframe.
#' @param na.action a function which indicates what should happen when the data contain 'NA's. Default is to fail.
#' @param alpha parameter of the one-step M-estimator. By default equals to 1.
#' @param cf consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
covW <- function (X, na.action = na.fail, alpha = 1, cf = 1 )
{
  X <- na.action(X)
  X <- as.matrix(X)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (p < 2)
    stop("'X' must be at least bivariate")
  Xmeans <- colMeans(X)
  di <- mahalanobis(X, Xmeans, cov(X))
  X.centered <- sweep(X, 2, Xmeans)
  v.tilde <- 1/n*cf * t( X.centered) %*% diag(di^alpha) %*%  X.centered

  return(v.tilde)
}

#' Mean Vector and Scatter Matrix Based on 4th Moments
#'
#' The function returns for some multivariate data the mean vector and the scatter matrix based on 4th moments.
#' @param x a numeric data matrix.
#'
#' @return
#' @export
#'
#' @examples
MeanCov4 <- function (x){
  list(location = colMeans(x), scatter = cov4(x))
}

#' Mean Vector and One-step M-estimator
#'
#' The function returns for some multivariate data the mean vector and the
#'  One-step M-estimator
#'
#'
#' @param x a numeric data matrix.
#' @param alpha parameter of the one-step M-estimator. By default equals to 1.
#' @param cf consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
MeanCovW <- function (x, alpha, cf){
  list(location = colMeans(x),
       scatter = covW(x, alpha = alpha, cf = cf))
}


#' Squared ICS Distances for Invariant Coordinates
#'
#' Computes the squared ICS distances, defined as
#' the Euclidian distances of the selected centered
#' components.
#'
#' @param dt_scores The dataframe of all Invariant Coordinates
#' @param index vector of integers indicating the indices of the components to select.
#'
#' @return
#' @export
#'
#' @examples
ics.distances.qr <- function (dt_scores, index = NULL)
{
  # if (class(object) != "ics2")
  #   stop("'object' must be of class 'ics2'")
  if (is.null(index))
    index <- 1:ncol(dt_scores)
  DIST <- rowSums(dt_scores[, index, drop = FALSE]^2)
  return(DIST)
}


#' Plot of ICS Distances
#'
#' @param df The initial nxp data frame
#' @param ind_outlier The vector of indexes of the outliers
#' @param k The indexes of the coordinates to keep
#' @param res_ICS_QR The output returned by ICS_QR functions
#' @param df_name The name of the data set
#' @param save Boolean TRUE or FALSE. By default FALSE if we do not want to save the plot.
#'
#' @return
#' @export
#'
#' @examples
plot.ics.distances <- function(df, res_ICS_QR, ind_outlier, k,
                               df_name = "ICSD", save = FALSE ){
  df_scores <- data.frame(ID = 1:nrow(df),
                          Defect = factor(ifelse(rownames(df)== ind_outlier, 1,0)),
                          Z = ics.distances.qr(res_ICS_QR$Z, index = k))
  # Create the y-axis label based on the combination of scatter matrices and
  # the number of selected components
  k_details <- paste(length(k),ifelse(length(k)==1, k,
                                      paste(range(k), collapse = ":")),
                     sep = ", ")
  y_axis <- bquote("ICSD"[ ~ k == .(k_details) ] ^2 )
  title <- title_alpha(alpha_val = res_ICS_QR$alpha)


  p <- df_scores %>%
    ggplot(aes(x = ID, y =  Z)) +
    geom_point(aes(color = Defect, shape = Defect), alpha = 0.7, size = 2.5) +
    scale_colour_manual(values = c("0" = "darkgrey", "1" = "red2")) +
   # scale_colour_colorblind()+
    theme_minimal() +
    labs(x = "Observation Number", y = y_axis, title = title) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20))

  if(save){
    save_pdf(title, p, df_name, width = 12)
  }else{
    return(list(p=p, title=title))
  }

}


#' Plot of ICS eigenvalues
#'
#' @param res_ICS_QR The output returned by ICQSR functions
#' @param df_name The name of the data set
#' @param save Boolean TRUE or FALSE. By default FALSE if we do not want to save the plot 
#'
#' @return
#' @export
#'
#' @examples
plot.ics.eigenvalues <- function(res_ICS_QR, df_name = "eigenvalues", save = FALSE){
  dt_eigenvalues <- data.frame(Index =  1:ncol(res_ICS_QR$Z),
                               Eigenvalues =  res_ICS_QR$eigenvalues)
  title <- title_alpha(alpha_val = res_ICS_QR$alpha)
  p <- dt_eigenvalues %>%
    ggplot( aes(x =Index, y =  Eigenvalues)) +
    geom_point() +
    theme_minimal() +
    labs( title = title) +
    theme(plot.title = element_text(hjust = 0.5))

  if(save){
    save_pdf(title, p, df_name)
  }else{
    return(p)
  }
}


#' Return the name of the combination of scatters based on alpha
#'
#' @param alpha_val The value of alpha
#'
#' @return
#' @export
#'
#' @examples
title_alpha <- function(alpha_val){
  if(alpha_val == 1){
    title = "Cov - Cov4"
  }else if (alpha_val == -1){
    title = "Cov - CovAxis"
  }else{
    title =   bquote("Cov - CovW"[ ~ alpha ~" "== " "~.(alpha_val) ] )
  }

  title
}


#' GGpairs save
#'
#' @param crabs crabs data only
#' @param res_ICS_QR The output returned by ICQSR functions
#' @param df_name The name of the data set
#' @param save Boolean TRUE or FALSE. By default FALSE if we do not want to save the plot
#'
#' @return
#' @export
#'
#' @examples
save_ggpairs <- function(crabs, res_ICS_QR, df_name, save = FALSE){
  # Pre process
  dt <- data.frame(sp = crabs$sp,
                   sex = crabs$sex,
                   group = paste(crabs$sp, crabs$sex, sep = "_"),
                   res_ICS_QR$Z)

  title <- title_alpha(alpha_val = res_ICS_QR$alpha)


  p <-  GGally::ggpairs(dt, aes(colour = group,
                                alpha = 0.4),
                        columns = 4:ncol(dt),
                        columnLabels = paste0("ICQR", 1:ncol(res_ICS_QR$Z)),
                        title =  title)
  if(save){
    save_pdf(title, p, df_name)
  }else{
    return(p)
  }

}

#' Save pdf
#'
#' @param title The combination of ICS scatters used.
#' @param p The plot to save.
#' @param df_name The name of the data set
#' @param width The width of the pdf. By default 7.
#'
#' @return
#' @export
#'
#' @examples
save_pdf <- function(title, p, df_name, width = 7){
  pdf_title <- janitor::make_clean_names(paste(title, collapse = ""),
                                         replace = c(`-0` = "neg_0"))
  pdf(paste0("figures/", df_name, "_", pdf_title, ".pdf"), width = width)
  print(p)
  dev.off()
}

#' Values of the mean vector to simulate a matrix
#' with k for condition number
#'
#' @param k The condition number
#' @param p The number of values to simulate
#'
#' @return A named vector with p mean values
#' @export
#'
#' @examples
k_sim <- function(k, p){
  if(k == 0){
    val <- rep(1,p)
  }else if (k=="ex"){
    val <- c(10^(-12), 10^(-3), 1, 10^6)
  }else{
    max_power <- round(k/2,0)
    min_power <- k - max_power
    all_val <- seq(max_power, -min_power)
    ind_val <- sample(2:(length(all_val)-1),p-2)
    # I want the minimum first
    val <- 10^all_val[c(length(all_val), ind_val, 1)]

  }

  names(val) <- letters[1:p]
  val
}

#' Simulations
#'
#' @param n Number of "normal" observations
#' @param n_outlier  Number of outliers
#' @param N_tot Total number of observations
#' @param p Number of variables
#' @param sigma  pxp Scatter matrix of "normal" observations
#' @param mean p-Vector of means of "normal" observations
#' @param sigma_outlier pxp Scatter matrix of outlying observations
#' @param mean_outlier p-Vector of means of outlying observations
#'
#' @return A N_totxp matrix. The first n rows are "normal" observations.
#' @export
#'
#' @examples
simulations = function(n = 980, n_outlier = 20,  p,
                       sigma = diag(c(1,rep(4,p-1))), mean =  rep(0,p),
                       sigma_outlier = diag(c(1,rep(4,p-1))), mean_outlier =  c(6,rep(0,p-1)) ){

  # We simulate the observations from the majority group
  x = rmvnorm(n = n, mean = mean, sigma = sigma)

  # We simulate the outliers
  x_outlier = rmvnorm(n = n_outlier, mean = mean_outlier, sigma = sigma_outlier)

  # We create our dataset with the first n_outlier observations as outliers
  X_ini = rbind( x, x_outlier)

  return(X_ini)

}

#' Simulation of ICA model
#'
#' @param n The number of observations to simulate
#'
#' @return
#' @export
#'
#' @examples
simulation_ICA = function(n){
  s1 <- rnorm(n)
  s2 <- rt(n,df=5) * sqrt(0.6)
  s3 <- runif(n, -sqrt(3),sqrt(3))
  s4 <- VGAM::rlaplace(n, location = 0, scale = 1/sqrt(2))

  S <- cbind(s1,s2,s3,s4)

   return(S)

}



#' Computing ICS eigenvalues through ICSEigen and ICSQR after
#' multiplying each column by different scales
#'
#' @param X_ini The initial nxp dataframe
#' @param k The vector of values to multiply each variable
#' @param alpha The value of alpha to compute CovWalpha
#' @param cf  The consistency factor of the one-step M-estimator. By default equals to 1
#'
#' @return
#' @export
#'
#' @examples
eigenvalues_all_methods <- function(X_ini, k, alpha, cf = 1){
  print(k)
  mean_all <- mean_all_save[,k]
  print(mean_all)
  X <-  sweep(X_ini, 2, mean_all, "*")
  print(colMeans(X))
  print(format(kappa(X, exact = TRUE), scientific = TRUE,
               digits = 2))

  # Theoretical --------------
  theoretical <- alpha==1 & "ICSTheory" %in% all_methods
  if(theoretical){
    gamma11 <- (1-delta)^2*(1-eps)*eps+1
    gamma21 <-((1-eps)*(3+6*(1-delta)^2*eps^2+
                          (1-delta)^4*eps^4)+
                 eps*(3+6*(1-eps)^2*(1-delta)^2+
                        (1-eps)^4*(1-delta)^4))

    ICS_th_eigenvalues <- c(as.vector((gamma21/(gamma11^2)+(p-1))/(p+2)),
                            rep(1,p-1))
    if(cf==1){ ICS_th_eigenvalues <- ICS_th_eigenvalues *(ncol(X)+2)}
    ICS_th_eigenvalues
  }

  # ICSEigen -------

  ics2_eigenvalues <- rep(NA,p)
  try( ics2_eigenvalues <-   ics2(X, S1 = MeanCov,
                                  S2 = MeanCovW, S2args = list(alpha = alpha, cf = cf)
  )@gKurt)


  # ICSQR --------------
  ICS_QR_eigenvalues <- rep(NA,2)
  try( ICS_QR_eigenvalues <- ICS_QR_full_rank(Xt = X, alpha = alpha, cf = cf)$eigenvalues)

  if(theoretical){
    return(list(ics2_eigenvalues,
                ICS_QR_eigenvalues, ICS_th_eigenvalues))
  }else{
    return(list(ics2_eigenvalues,
                ICS_QR_eigenvalues))
  }

}




#' Plot of comparison of eigenvalues based on different computations
#'
#' @param all_methods The names of the methods used to compute the eigenvalues
#' @param eigenvalues_all_dt The list of all eigenvalues obtained with the different methods
#' with the eigenvalues_all_methods function.
#' @param alpha  parameter of the one-step M-estimator. By default equals to 1.
#' @param df_name The name of the data set
#' @param save Boolean TRUE or FALSE. By default FALSE if we do not want to save the plot
#' @param df_name The name of the data set
#' @param nb The number of initial variables.
#' @param k_all The vector of condition numbers
#'
#' @return
#' @export
#'
#' @examples
dt_eigenvalues_plot <- function(all_methods, eigenvalues_all_dt, nb,
                                alpha, df_name = "", save = FALSE, k_all){
  
  theoretical <- alpha==1 & "ICSTheory" %in% all_methods
  if(theoretical){
    ics_theory <-  eigenvalues_all_dt[[1]][[3]]
  }
  all_methods <- all_methods[ all_methods != "ICSTheory"]
  dt_all_list <- lapply(1:length(all_methods), function(j){
    dt <- data.frame(t(sapply(1:length(eigenvalues_all_dt),
                              function(i) eigenvalues_all_dt[[i]][[j]])))
    dt <- cbind(methods = all_methods[j],
                K = 0:(nrow(dt)-1),
                dt)
    dt
  })
  dt_all <- data.table::rbindlist( dt_all_list, use.names = TRUE,
                                   fill = TRUE)


  dt2 <- reshape(dt_all, varying = list(3:(nb+2)), v.names = "value",
                 timevar = "eigenvalue", direction = "long")
  title <- title_alpha(alpha_val = alpha )
  p <- dt2 %>% ggplot(aes(x = K, y = value,
                         color = factor(eigenvalue))) +
    {if(theoretical) geom_hline(yintercept = ics_theory,
                                linetype = "dashed",
                                color = "darkgrey") }+
    geom_point() +
    geom_line() +

    ylim(0,NA)+
    theme_minimal() +
    labs(  title = title,
           color = "Eigenvalues",
           x = "Condition number (in log10)"
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 20),
          legend.position="top") +
    scale_colour_colorblind() +
    facet_grid(~methods)
p

  if(save){
    save_pdf(title, p, df_name, width = 7*length(all_methods))
  }else{
    return(p)
  }


}

