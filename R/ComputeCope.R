#' Compute CoPE sets.
#' 
#' Computes CoPE sets for the data Y using the algorithm from Sommerfeld, Sain
#' and Schwartzman (2015).
#'
#' @param Z A list with components "x", "y" and "z". Here, x and y are the 
#'          coordinates of the grid on which the data is observed and z is an 
#'          array with dimensions c(length(x),length(y),n), containing the
#'          data. n is the number of observations.
#' @param level The level of interest.
#' @param X The design matrix of the linear model. If NULL, it is set to 
#'          matrix(rep(1,dim(Y)[3]),ncol=1) corresponding to i.i.d. data.
#' @param w Vector of weights. The target function is w^T * beta where beta is 
#'          the true parameter function. By default the target function is the
#'          first entry of beta.
#' @param correlation Type of correlation assumed for the spatially indexed 
#'                    indexed linear models. This is a string that is passed to
#'                    the function gls from the nlme package. Defaults to NULL 
#'                    which corresponds to i.i.d. errors.
#' @param corpar A list of parameters passed to the correlation function.
#' @param alpha The significance level. Inclusion for CoPE sets holds with 
#'              probability 1-alpha.
#' @param N The number of bootstrap realizations to generate for determining
#'          the threshold.
#' @param mu The true parameter function. Use the default NULL if unknown. 
#' @param mask Pixels outside the mask (i.e. where mask is ==NA) are ignored.
#'   
#' @return An object of class cope. A list containing the following
#' \itemize{
#'  \item{x, y: }{The grid coordinates from the input.}
#'  \item{mu, level, tau, X, N, alpha, mask: }
#'  {The corresponding values from the input.}
#'  \item{mu_hat, norm_est: }{The estimatot for mu and its normalized version.}
#'  \item{a_MB, a_MB_true, a_Tay, a_Tay_true: }{Thresholds for the CoPE sets 
#'  determined using the multiplier bootstrap or Taylor's method and the 
#'  estimated or the true contour, respectively.}
#'  \item{incl_MB, incl_MB_true, incl_Tay, incl_Tay_true: }{Booleans indicating 
#'  whether inclusion of the excursion set in the upper CoPE set and inclusion
#'  of the lower CoPE set in the excursion set holds, when CoPE sets are
#'  determined by a_MB, a_MB_true, a_Tay or a_Tay_true, respectively. Only 
#'  available if mu is given.}
#' }
#' 
#' @export
#' 
#' @references M. Sommerfeld, S. Sain and A. Schwartzman. Confidence regions for 
#'             excursion sets in asymptotically Gaussian
#'             random fields, with an application to climate. Preprint, 2015. 
#'        
#' @examples
#' # An example using the ToyNoise and ToySignal of this package.
#' n = 30
#' Data = ToyNoise1(n = n)
#' Data$z = Data$z + rep(ToySignal()$z, n)
#' CopeSet = ComputeCope(Data,level=4/3, mu=ToySignal()$z)
#' PlotCope(CopeSet)
ComputeCope = function(Z,level,
                X=NULL,
                w = NULL,
                correlation = NULL,
                corpar = NULL,
                groups = NULL,
                alpha=0.1,
                N=1000,
                mu=NULL,
                mask=NULL){
  x = Z$x
  y = Z$y
  Y = Z$z
  n = dim(Y)[3]
  nloc <- length(x) * length(y)
  
  if(is.null(X)){
    X <- matrix(1, n, 1)
    w <- matrix(1, 1, 1)
  } 
  p <- ncol(X)
  
  if(is.null(groups)){
    groups <- rep(1, n)
  }
  
  
  # Compute the inverse of the square-root of a positive definite matrix.
  invsqrtm <- function(A){
    E <- eigen(A)
    U <- E$vectors
    D <- diag(E$values)
    U %*%  diag(1 / sqrt(E$values)) %*% t(U)
  }
 
  # For each location compute the GLS estimate, the de-correlated residuals and 
  # vabs and the normalized estimator.
  deR <- array(0, c(length(x), length(y), n)) # De-correlated residuals.
  vabs <- matrix(0, length(x), length(y)) # Absolute value of v.
  norm_est <- matrix(0, length(x), length(y)) # Normalized estimator.
  mu_hat <- matrix(0, length(x), length(y)) # Estimate of the function of 
                                            # interest.
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      ytemp <- Y[i, j, ]
      df <- data.frame(cbind(ytemp = ytemp, X, groups = groups))
      df <- df[order(groups), ]
      groups <- sort(groups)
      
      fo <- paste(names(df)[1], "~", 
                  paste(names(df)[-c(1, p + 2)], collapse=" + "), "-1")
      
      if(is.null(correlation)){
       model <- nlme::gls(formula(fo), data = df)       
      }else{
       model <- nlme::gls(formula(fo), data = df,  
                          correlation = do.call(get(correlation),  
                                                c(corpar, form = ~1|groups)))       
      }

      mu_hat[i, j] <- t(w) %*% model$coefficients
      
      if(!is.null(correlation)){
        cM <- corMatrix(model$modelStruct$corStruct)
        if(!is.list(cM)) cM <- list(cM)
        invsqrtmOmega <- Matrix::as.matrix(Matrix::bdiag(lapply(cM, invsqrtm)))
        deR[i, j, ] <- invsqrtmOmega %*% model$residuals 
        deR[i, j, ] <- deR[i, j, ] / sd(deR[i, j, ])
      } else{
        deR[i, j, ] <- model$residuals 
        deR[i, j, ] <- deR[i, j, ] / sd(deR[i, j, ])
      }
      
      vabs[i, j] <- sqrt(t(w) %*% model$varBeta %*% w)
      
      norm_est[i, j] <- (mu_hat[i, j] - level) / vabs[i, j] 
    }
    print(round(i / length(x) * 100))
  }
  
  #Ignore values not on mask.
  if(!is.null(mask)){
    mu_hat[is.na(mask)] = level - 10
    norm_est[is.na(mask)] <- NA
   } 
  
  
  #Compute a using Multiplier Bootstrap.
  a_MB = quantile(MBContour(x=x, y=y, R=deR, N=N, f=mu_hat, level=level),
                  probs=1-alpha,type=8)
  #Compute a using Taylors method.
  Tay_fun = TaylorContour(x=x, y=y, f=mu_hat, level=level, R=deR)
  a_Tay = 0
  while(2*Tay_fun(a_Tay)>alpha) a_Tay = a_Tay + 0.01  #The two is for the absolute value.
  
  #Compute threshold a with the true boundary if available.
  if(is.null(mu)) {a_MB_true = NA; a_Tay_true = NA} else{
    #Compute a using Multiplier Bootstrap.
    a_MB_true = quantile(MBContour(x=x,y=y,R=deR,N=N,f=mu,level=level),
                         probs=1-alpha,type=8)
    #Compute a using Taylors method.
    Tay_fun = TaylorContour(x=x,y=y,f=mu,level=level,R=deR) #The two is for the absolute value.
    a_Tay_true = 0
    while(2*Tay_fun(a_Tay_true)>alpha) a_Tay_true = a_Tay_true + 0.01
  }
  
  
  #Determine whether inclusion holds.
  if(is.null(mu)){incl_MB = NA; incl_Tay = NA; incl_MB_true = NA; incl_Tay_true = NA} else{
    incl_MB = SubsetContour(x,y,norm_est,a_MB,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_MB)
    incl_Tay = SubsetContour(x,y,norm_est,a_Tay,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_Tay)
    incl_MB_true = SubsetContour(x,y,norm_est,a_MB_true,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_MB_true)
    incl_Tay_true = SubsetContour(x,y,norm_est,a_Tay_true,mu,level) & SubsetContour(x,y,mu,level,norm_est,-a_Tay_true)
  }
  
  result = structure(
    list(x=x,y=y,mu=mu,level=level,mu_hat=mu_hat,w = w,X=X,alpha=alpha,N=N,norm_est=norm_est, 
         a_MB=a_MB,a_Tay=a_Tay,a_MB_true=a_MB_true,a_Tay_true=a_Tay_true,
         incl_MB=incl_MB, incl_Tay=incl_Tay,incl_MB_true=incl_MB_true,incl_Tay_true=incl_Tay_true, mask=mask), class = "cope")
}

#' Plots CoPE sets. 
#'
#' @param cope An object of class cope to be plotted.
#' @param plot.taylor  Boolean indicating whether the CoPE sets with the threshold 
#'                      obtained by Taylor's method should be plotted. Default is
#'                      FALSE. 
#' @param use.true.function  Boolean indicating whether the threshold obtained 
#'                            from the true function should be used. Default is 
#'                            FALSE.  
#' @param map If TRUE plot the cope set on a map of the world. The coordinates 
#'            in this case are interpreted as longitude and latitude.  
#' @param ... Additional graphical parameters passed to fields::image.plot.
#' @export
PlotCope = function(cope,plot.taylor=FALSE, use.true.function = FALSE, map=FALSE, ...){
  
  if(map){
    plot.function = ImageMap
  } else{
    plot.function = fields::image.plot
  }
  
  if(use.true.function & is.null(cope$mu)){
    print("True function is not available. Estimated boundary will be used.")
  } 
  
  x = cope$x
  y = cope$y
  
  if(use.true.function & !is.null(cope$mu)){
    if(plot.taylor) a = cope$a_Tay_true else a = cope$a_MB_true
    plot.function(x,y,cope$mu_hat,horizontal=FALSE, mask=cope$mask, ...)
    DrawContour(x,y,cope$mu,level=cope$level,col="purple")
    DrawContour(x,y,cope$norm_est,level=a,col="darkred")
    DrawContour(x,y,cope$norm_est,level=-a,col="darkgreen")
  } else{
    if(plot.taylor) a = cope$a_Tay else a = cope$a_MB
    plot.function(x,y,cope$mu_hat,horizontal=FALSE, mask=cope$mask, ...)
    if(!is.null(cope$mu)){
      DrawContour(x,y,cope$mu,level=cope$level,col="purple")
    } else{
      DrawContour(x,y,cope$mu_hat,level=cope$level,col="purple")
    }
    DrawContour(x,y,cope$norm_est,level=a,col="darkred")
    DrawContour(x,y,cope$norm_est,level=-a,col="darkgreen")
  }
}