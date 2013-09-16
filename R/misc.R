LL.alt<-function(mat,model,phenotype.dat)
{
  if(ncol(mat)==8){
    qtlmodel<-as.formula(paste(deparse(model),"+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]",sep=""))
  }else{
    if(ncol(mat==16))
    {
      qtlmodel<-as.formula(paste(deparse(model),
                                 "+mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]
                                 +mat[,9]+mat[,10]+mat[,11]+mat[,12]+mat[,13]+mat[,14]+mat[,15]
                                 ",sep=""))
      
    }else{
      stop("problem with creating genotype matrix")  
    }#catch else
  }#initial else
  logLik(lm(qtlmodel,data=phenotype.dat))/log(10)
}


LL.multi<-function(mat,model,pheno)
{
  if(ncol(mat)==8){
    qtlmodel<-as.formula("pheno~mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]")
  }else{
    if(ncol(mat==16))
    {
      qtlmodel<-as.formula("pheno~mat[,1]+mat[,2]+mat[,3]+mat[,4]+mat[,5]+mat[,6]+mat[,7]
                           +mat[,9]+mat[,10]+mat[,11]+mat[,12]+mat[,13]+mat[,14]+mat[,15]")
      
    }else{
      stop("problem with creating genotype matrix")	
    }#catch else
  }#initial else
  logLik.multi(lm(qtlmodel))/log(10)
  
}

logLik.multi <- function(object, REML = FALSE, ...)
{
  
  all.val<-numeric(length=0)
  for(j in 1:ncol(object$residuals))
  {
    res <- object$residuals[,j] # not resid(object) because of NA methods
    p <- object$rank
    N <- length(res) 
    if(is.null(w <- object$weights)) {
      w <- rep.int(1, N)
    } else {
      ## this is OK as both resids and weights are for the cases used
      excl <- w == 0    	# eliminating zero weights
      if (any(excl)) {
        res <- res[!excl]
        N <- length(res)
        w <- w[!excl]
      }
    }
    N0 <- N
    if(REML) N <- N - p
    val <- .5* (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) +
                                     log(sum(w*res^2))))
    if(REML) val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
    attr(val, "nall") <- N0 # NB, still omits zero weights
    attr(val, "nobs") <- N
    attr(val, "df") <- p + 1
    class(val) <- "logLik"
    all.val<-rbind(all.val,val)
  }#for close
  all.val
}
