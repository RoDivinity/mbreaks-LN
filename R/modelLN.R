
#create new S3 class (building block of everything):
#S3 class to capture the structural break model, will have the following
#properties:
#   numbr: number of breaks (not maximum, but a specify number of breaks); other S3 class can 
#inherit a different meaning of numbr later which is different, the maximum number of break you want
# numbr is numeric
#   mform: special type of class formula, which can understand the syntax.
# only contains breaks in parameters, not variance
# `y~z|x`, where x could be empty



#The model should be:
#Input is a formula, a data


#initialization of class sbm (only care about estimation of structural break model given
# number of breaks, not concerning with break selection. Break selection model should be
#subclass inherits many of the methods below from this block model/class)
#this class will not use trimming (no relationship to testing?) -> to simplify, this one
# might not have testing functions available for number of break dates. to test for 
# number of breaks, we will have second class that is capable of this
# under control, we have: 
# robust; prewhit; hetdat; hetvar; hetomega; hetq?
#' @export
sbm <- function(mform, data, m=2, h=10, control=NULL,  ...){
  
  envForm <- attr(mform,'.Environment')
  #check validity of formula in structural break model
  sbm.formula = sbm.check_formula(mform)
  
  
  #check for constant/intercept in the model + check for duplicates regressors
  vars = check_const(sbm.formula$zvars,sbm.formula$xvars)
  zvars = vars$zt
  xvars = vars$xt
 
  #add special case, where no real xvars is specified but users put only -1.
  if ( '-1' %in% xvars){
    xvars <- xvars[-which(xvars == '-1')]
  }
  
  #default list of controls/options
  con = list('robust' = TRUE,'hetomega'=TRUE, 'hetdat' = TRUE, 'hetq'=TRUE,'hetvar'=TRUE,'prewhit'=FALSE)

  #check for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  
  
  
  ## process and transform data to appropriate y,z,x matrices
  #fix option = na.fail for handling missing values for now
  #since this is time series, missing observations cannot be tolerant
  if (missing(data)){
  sbm.frame = model.frame(sbm.formula$tform, na.action = na.fail)}
  else{
    sbm.frame = data
  }

  #match data frame with variables  
  y_ind = match(sbm.formula$yvar,colnames(sbm.frame))

  
  y = as.matrix(sbm.frame[,y_ind,drop=FALSE])
  bigt = dim(y)[1]
  #check if any observations
  if (bigt == 0) {stop('No observations in the model')}
  #check if any z regressors
  if (length(zvars)==0){ stop('No z regressors. Use `lm` instead')}
  z = matrix(0,0,0)
  x = matrix(0,0,0)
  if (attr(zvars,'const') == 1){
    z = matrix(1,bigt,1)
    colnames(z)[1] <- '(Intercept)'
  }
  if (attr(xvars, 'const') == 1){
    x = matrix(1,bigt,1)
    colnames(x)[1] <- '(Intercept)'
  }
  
  #add additional regressors from frame
  z_ind = match(zvars,colnames(sbm.frame[,-1]))
  x_ind = match(xvars,colnames(sbm.frame[,-1]))
  #form regressor matrices (x is optional)
  z=cbind(z,as.matrix(sbm.frame[,z_ind[which(!is.na(z_ind))],drop=FALSE]))
  if (identical(x,matrix(0,0,0)) && identical(x_ind, integer(0))){}else{
  x=cbind(x,as.matrix(sbm.frame[,x_ind,drop=FALSE]))}
  
  
  #compute other parameters to fit in arguments of aux functions
  p = dim(x)[2]
  q = dim(z)[2]
  
  #check validity of minimum segment length options &
  #check validity of number of breaks
  #check_segOpt(h,m,bigt) 
  
  #estimate the break date with given number of break
  if (p==0){
    datevec <- dating_purescSSR(y, z, m, h)$datevec 
  }else{
    datevec <- dating_partscSSR(y, z, x, m, h)$datevec 
  }
  est_m = m
  brdate <- as.matrix(datevec[,est_m])
  #estimate the sbm model parameters with estimated break date (for inference)
  sbm.est <- estim(y,z,x,est_m,brdate,con$robust,con$prewhit,con$hetomega,con$hetq,con$hetdat,con$hetvar)
  out <- c()
  #general description of the model
  rhs_z_ref <- paste(vars$zt,collapse='+')
  if (p == 0){rhs_ref <- rhs_z_ref}else{
    rhs_x_ref <- paste(vars$xt,collapse='+')
    rhs_ref <- paste(rhs_z_ref,rhs_x_ref,sep='|')}
  form_ref <- paste(sbm.formula$yvar,rhs_ref,sep='~')
  
  #group data into one separate field model.sbm
  model.sbm <- c()
  model.sbm$y <- y
  model.sbm$x <- x
  model.sbm$z <- z
  attr(model.sbm, 'dep.var') <- colnames(y)
  attr(model.sbm, 'xterms') <- colnames(x)
  attr(model.sbm, 'zterms') <- colnames(z)
  #group description of the model into field description
  description <- c()
  out.formula <- formula(form_ref)
  attr(out.formula,'.Environment') <- envForm
  description$formula <- form_ref
  attr(description, 'formula') <- form_ref
  attr(description, 'MAXbreaks') <- m
  attr(description, 'min.segLength') <- h
  if (p==0){attr(description, 'type') <- 'Pure'}
  else{attr(description, 'type') <- 'Partial'}
 
 
  # break date estimates and CI stored as bound (for inference)
  breakdate <- c()
  breakdate$date <- brdate
  breakdate$bound <- sbm.est$bound
  attr(breakdate, 'numbreaks') <- est_m #this denotes differently
  #to reflect it is estimated by IC/test (for later more complex class)
  #could be more way of estimating confidence interval for break date later
  attr(breakdate, 'procedure') <- 'Asymptotic'
  # regressors estimates and CI stored as stdev (for inference)
  coefficients <- c()
  sbm.est$coef <- as.vector(sbm.est$coef)
  coefficients$coef.z <- sbm.est$coef[1:((est_m+1)*q)]
  coefficients$SE.cor.z <- sbm.est$SE_correct[1:((est_m+1)*q)]
  if (p == 0){
    coefficients$coef.x <- NULL
    coefficients$SE.cor.x <- NULL}
  else{
    coefficients$coef.x <- sbm.est$coef[((est_m+1)*q+1):((est_m+1)*q+p)]
    coefficients$SE.cor.x <- sbm.est$SE_correct[((est_m+1)*q+1):((est_m+1)*q+p)]}
  
  
  attr(coefficients, 'xregs') <- colnames(x)
  attr(coefficients, 'zregs') <- colnames(z)
  attr(coefficients, 'numx') <- p
  attr(coefficients, 'numz') <- q
  #add backs all above fields
  out$description <- description
  out$model <- model.sbm
  out$breakdate <- breakdate
  out$coefficients <- coefficients
  class(out) <- 'sbm'
  #post-regression 
  out$deviance = sbm.est$SSR
  out$residuals = sbm.est$resid
  out$fitted.values = sbm.est$fitted
  # keep tracks of options used
  out$controls = con
  #format the results
  out$table <- compile(out)
  return(out)
}


#Lists of methods to do
#1) summary
#2) print (format as suggested)
#3) plot*** (important)
#4) coef (must be able to split out estimated coefs on z and x (optional))
#5) brdate (get the estimated break date out, in vector and other form for display)
#6) model.frame (need to rewrite this), which extracting model frame from a formula/fit
# (Unnecessary, can reuse base R)
#7) resid (return already computed residuals of the model)


#check object class
#' @export 
is.sbm <- function(x) inherits(x, 'sbm')

##generic methods declaration


#lists of generic functions for class sbm
#' @export summary
summary <- function(x, ...){
  UseMethod('summary',x)
}
#' @export print
print <- function(x, ...){
  UseMethod('print',x)
}
#' @export plot
plot <- function(x, ...){
  UseMethod('plot',x)
}
#' @export coef
coef <- function(x, ...){
  UseMethod('coef',x)
}
#' @export residuals
residuals <- function(x,...){
  UseMethod('residuals',x)
}
#' @export fitted
fitted <- function(x,...){
  UseMethod('fitted',x)
}

#' @export
brdate <- function(x,...){
  UseMethod('brdate',x)
}

#' @export compile
compile <- function(x,...){
  UseMethod('compile',x)
}

#' @export getOption
getOption <- function(x,...){
  UseMethod('getOption',x)
}


#specific S3 functions for class sbm
#' @export
compile.sbm <- function(x, digits = 3,...){
    ## format date estimation
    bound95 = c()
    bound90 = c()
    coln = c()
    numbreaks = attr(x$breakdate, 'numbreaks')
    for (i in 1:numbreaks){
      coln = c(coln, paste('Break',i,sep=''))
      bound95 = c(bound95,paste('(',x$breakdate$bound[i,1],',',x$breakdate$bound[i,2],')',sep=''))
      bound90 = c(bound90,paste('(',x$breakdate$bound[i,3],',',x$breakdate$bound[i,4],')',sep=''))
    }
    date_tab = data.frame(date = t(x$breakdate$date),stringsAsFactors = FALSE)
    date_tab = rbind(date_tab,bound95)
    date_tab = rbind(date_tab,bound90)
    colnames(date_tab) = coln
    rownames(date_tab) = c('Date','95% CI','90% CI')
    
    ## format z regressors coefficients (varying estimates over regimes)
    q = attr(x$coefficients, 'numz')
    rnameRS = attr(x$coefficients, 'zregs')
    cnameRS = c()
    coefRS = c()
    for (i in 1:(numbreaks+1)){
      cnameRS = cbind(cnameRS,paste('Regime',i))
    }
    for (j in 1:q){
      coefRSj = c()
      for (i in 1:(numbreaks+1)){
        coef_val = format(round(x$coefficients$coef.z[(j-1)*(numbreaks+1)+i],digits),nsmall=digits)
        coef_std = paste('(',format(round(x$coefficients$SE.cor.z[(j-1)*(numbreaks+1)+i],digits),nsmall=digits),')',sep='')
        coefRSj = cbind(coefRSj,paste(coef_val,coef_std,sep = ' '))
      }
      coefRS = rbind(coefRS,coefRSj)
    }
    
  rnameRSf = paste(rnameRS,'(SE)',sep=' ')
  coef_tabRS=data.frame(coef = coefRS)
  rownames(coef_tabRS) = rnameRSf
  colnames(coef_tabRS) = cnameRS
  
  ## format x regressors coefficients (constant estimates over regimes)
  p = attr(x$coefficients, 'numx')
  if(p == 0){coef_tabRW = NULL}
  else{
    rnameRW = attr(x$coefficients, 'xregs')
    cnameRW = 'All regimes'
    coefRW = c()
    for (j in 1:p){
      coef_val = format(round(x$coefficients$coef.x[j],digits),nsmall=digits)
      coef_std = paste('(',format(round(x$coefficients$SE.cor.x[j],digits),nsmall=digits),')',sep='')
      coefRW = cbind(coefRW,paste(coef_val,coef_std,sep=' '))}
    
    rnameRWf = paste(rnameRW,'(SE)',sep=' ')
    coef_tabRW=data.frame(coef = coefRW)
    colnames(coef_tabRW) = rnameRWf
    rownames(coef_tabRW) = cnameRW}
  
  table = list('date_tab' = date_tab,'zregs_tab' = coef_tabRS, 'xregs_tab' = coef_tabRW)
  return(table)
}

#' @export
# why this one do not appear when call S3 object to console?
print.sbm <- function(x, ...){
  cat(paste(gettextf('\n%s structural break model with %d breaks',
              attr(x$description,'type'),attr(x$breakdate,'numbreaks')),'\n'))
  cat(paste('Model specification: ',attr(x$description,'formula'),"\n"))
  if (attr(x$breakdate,'numbreaks') < 1) {
   cat(paste('\nNo break estimated'))
  }else{
   print(x$table$date_tab,quote=FALSE)}
  cat("\n")
}


#' @export
summary.sbm <- function(x, ... ){
  cat(paste(gettextf('%s structural break model with %d estimated breaks',
               attr(x$description,'type'),attr(x$breakdate,'numbreaks')),'\n'))
  
  cat(paste('Model specification: ',attr(x$description,'formula'),"\n"))
  
  cat('\nEstimated date:\n')
  print(x$table$date_tab,quote=FALSE)
  
  cat('\nEstimated regime-specific coefficients:\n')
  print(x$table$zregs_tab,quote=FALSE)
  
  if(attr(x$coefficients, 'numx') == 0) {cat('\nNo full sample regressors\n')}
  else{
    cat('\nEstimated full-sample coefficients:\n\n')
    print(x$table$xregs_tab,quote='FALSE')}
  cat('\nMinimum SSR =',
    format(round(x$deviance,3),nsmall=3),'\n')
  invisible(x)
}


#' @export
brdate.sbm <- function(x,...) x$breakdate$date

#' @export
coef.sbm <- function(x,...){ 
  out <- c()
  numbreaks = attr(x$breakdate, 'numbreaks')
  p = attr(x$coefficients, 'numx')
  q = attr(x$coefficients, 'numz')
  cname.z = c()
  coef.z = matrix(0L,q,numbreaks+1)
  rname.z= attr(x$coefficients, 'zregs')
  for (i in 1:(numbreaks+1)){
    cname.z = cbind(cname.z,paste('Regime',i))
  }
  print(x$coefficients$coef.z)
  for (j in 1:q){
      coef.z[j,] = x$coefficients$coef.z[((j-1)*(numbreaks+1)+1):(j*(numbreaks+1))]}
  colnames(coef.z) <- cname.z
  rownames(coef.z) <- rname.z
  coef.x = x$coefficients$coef.x
  names(coef.x) <- attr(x$coefficients, 'xregs')
  out$coefficients.z <- coef.z
  out$coefficients.x <- coef.x
  out
}

#' @importFrom graphics lines plot legend segments abline
#' @export
plot.sbm <- function(x,caption = NULL,xlab = NULL,ylab = NULL,bound=0.95,null.model = TRUE,start = NULL,...){
  est_m <- attr(x$breakdate, 'numbreaks')
  dev.hold()
  #plot x-axis based on provided start date (need to implement later).
  if(!is.null(start)) {}
  #comparison between structural break vs no break
  y <- unclass(x$model$y)
  bigt <- dim(y)[1]
  fity <- x$fitted.values
  x_t <- seq(1,dim(y)[1],1)
  range_y <- max(y)-min(y)
  if(is.null(ylab)){ylab <- colnames(y)}
  if(is.null(xlab)){xlab <- 'Time'}
  graphics::plot(x_t,y,type='l',col="black", xlab=xlab,ylab=ylab, 
       ylim=c(min(y)-range_y/10,max(y)+range_y/10),lty=1)
  #plot fitted values series for break model
  graphics::lines(x_t, fity,type='l', col="blue",lty=2)
  
  if(null.model){ #turn on plot for no break model
    zreg = x$model$z
    if(attr(x$coefficients,'numx') == 0){xreg=matrix(0,bigt,0)}else{
    xreg = x$model$x}
    fixreg = cbind(xreg,zreg)
    fixbeta = olsqr(y,fixreg)
    fity_fix = fixreg%*%fixbeta
    #plot fitted values series for null model
    graphics::lines(x_t, fity_fix,type='l', col="dark red",lty=2)
  }
  
  if(null.model){
    graphics::legend('topleft',legend=c("true",
                                        paste(est_m,'break(s)'),"0 break"),
           lty=c(1,2,2), col=c("black","blue","red"), ncol=1,bty='n')
  }else{
    #0,max(y)+range_y/10
    graphics::legend('topleft',legend=c("true",paste(est_m,'break(s)')),
           lty=c(1,2), col=c("black","blue"), ncol=1,bty='n')
  }
  #plot estimated dates + CIs based on bound option
  
  bigT <-length(y)
  for (i in 1:est_m){
    graphics::abline(v=x$breakdate$date[i,1],lty=2)}
  if (is.null(bound)){}
  else if (bound == 0.95){
    for (i in 1:est_m){
      lb_i = x$breakdate$bound[i,1]
      ub_i = x$breakdate$bound[i,2]
      if (lb_i < 0){lb_i = 0}
      if(ub_i>bigT){ ub_i=bigT}
      graphics::segments(lb_i,min(y)*(12+i/est_m)/10,ub_i,min(y)*(12+i/est_m)/10,lty=1,col='red')
    }
  }
  else if (bound == 0.9){
  for (i in 1:est_m){
    lb_i = x$breakdate$bound[i,3]
    ub_i = x$breakdate$bound[i,4]
    if (lb_i < 0){lb_i = 0}
    if(ub_i>bigT){ ub_i=bigT}
    graphics::segments(lb_i,min(y)*(12+i/est_m)/10,ub_i,min(y)*(12+i/est_m)/10,lty=1,col='red')
  }
  }
  dev.flush()
  invisible()
}

#' @export
#return list of options used on segments and variance structure
getOption.sbm <- function(x,...){
  #default list of controls/options
  con = x$controls
  cat(gettextf('robust = %s; hetdat = %s; hetvar = %s; hetomega = %s; hetq = %s',
           con$robust,con$hetdat,con$hetvar,con$hetomega,con$hetq))
}

#' @export
residuals.sbm <- function(x,...) x$residuals
#' @export
fitted.sbm <- function(x,...) x$fitted.values




###### checking validity of formula

## process formula when specifying mbreaks model
# Input: mform: a formula specific to mbreak models. Require lhs, and rhs includes
# z regressors (coefficients allowed to change across regimes, including constant) 
# and x regressors (coefficients are not changed across regimes, excluding constant)
# Output: a list containing: name of dependent variable y
#
# Need to check no overlapping x and z regressors (separate help function for easier
# maintainance)
# To specify constant mean across regimes, must explicitly specify 1 in x regressors
#         name of z regressors
#         name of x regressors (if any)
#         a transformed formula in basic linear regression later used to construct data frame
sbm.check_formula <- function(mform) {
 #expression for evaluation must be a formula
 if(!inherits(mform,"formula")){stop('The expression for model is not a formula.')}
 
 #must be a 3 terms formula, y ~ x | z; if not, automatically stop the model due to 
 #missing dependent variable y
 if (length(mform) == 2){ stop('Missing dependent variable y')}
 rhs <- if(length(mform) > 2) mform[[3]] 
 lhs <- if(length(mform) > 2) mform[[2]] 
  
 if(is.null(rhs)) {stop('There is no z and x regressors')}
  
 #check for correct syntax in RHS variables (only either 0 or 1 separator | is allowed)
 yvar <- lhs
 zvars <- list()
 xvars <- list()
 flagNoX <- FALSE
 if(length(rhs) < 3 || rhs[[1]] != '|'){
  flagNoX <- TRUE
  zvars_p <- rhs
  zvars <- get_terms(zvars_p)
  message('No x regressors in formula')} #the formula contains only z regrssors
 else {
  zvars_p = rhs[[2]]
  xvars_p = rhs[[3]]
  if(length(zvars_p)>1 && zvars_p[[1]] == '|'){stop('Too many separator `|` in RHS regressors. You can only have x and z regressors')}
  if(length(xvars_p)>1 && xvars_p[[1]] == '|'){stop('Too many separator `|` in RHS regressors. You can only have x and z regressors')}
  #zvars <- extract_terms(zvars_p)
  #xvars <- extract_terms(xvars_p)
  zvars <-get_terms(zvars_p)
  xvars <-get_terms(xvars_p)
  }
 
 #reconstruct to a linear regression formula to create corresponding data frame later
 
 if (flagNoX){
  mform_t <- . ~ .
  mform_t[[2]] <- yvar
  mform_t[[3]] <- zvars_p
 }
 else{
  mform_t <- . ~ . + .
  mform_t [[3]][[2]] <- zvars_p
  mform_t [[3]][[3]] <- xvars_p
  mform_t[[2]] <- yvar
 }
 return(list('yvar'=as.list(yvar),'zvars' = zvars,'xvars' = xvars,'tform' = mform_t))
}


check_overlap <-function(x,z){
  if (is.null(x)){
  }
  else{
    #get the type regressors with more variables
    xlen = length(x)
    zlen = length(z)
    if (xlen > zlen){
      lvars = x
      svars = z}
    else{
      lvars = z
      svars = x}
    if (any(as.vector(lvars)%in%as.vector(svars))){
     ov.flag = as.vector(lvars)%in%as.vector(svars)
     ov.vars = lvars[which(ov.flag==TRUE)]
     stop(gettextf('Variables %s are both z and x type regressors',
                   toString(ov.vars)))}
  }
}

#check if regression model has constant terms in z or x or none 
#+check for overlap
check_const <- function(zterms,xterms){
    zterms <- as.vector(zterms)
    xterns <- as.vector(xterms)
    if ('-1' %in% zterms){
      #if no intercept is specified in z
      #remove '-1' in z terms
      zterms <- zterms[-which(zterms=='-1')]
      #set correct property
      attr(zterms,'const') <- 0
      #check if no intercept is specified in x, otherwise constant assume no change
      if('-1' %in% xterms){
        xterms <- xterms[-which(xterms=='-1')]
        attr(xterms,'const') <- 0}
      else{attr(xterms,'const') <- 1}}
    else {
      #constant is prioritize in z regressors
      attr(zterms,'const') <- 1
      attr(xterms,'const') <- 0 
    }
    #check for overlapping terms between x and z regressors
    check_overlap(xterms,zterms)
   list(zt = zterms,xt = xterms)
}

#collect variables name in a term
get_terms <- function(terms){
  nt <- c()
  if (length(terms) == 1) {
    nt = terms
    class(nt) = 'list'}
  else{
    while(length(terms)>1 && terms[[1]] == '+'){
      nt <- c(terms[[3]],nt)
      terms <- terms[[2]]
    }
    #add last term on the right
    nt <- c(terms, nt)
  }
  nt <- unique(nt)
  nt
}
# 
# #test case
# testf = y~x1+x2+x3+x4
# testf_1sep = y~x1+x2|x3+x4
# testf_2sep = y~x1|x2|x3+x4
# testf_3sep = y~x1|x2|x3|x4
# testf_dsep = y~x1&x2|x3+x4
# testf_dsep2 = y~x1|x2+x3&x4
# testf_dsep3 = y~x1|x2|x3&x4
# 
# # #debug code cont
# x1 = runif(10)
# x2 = runif(10)
# x3 = runif(10)
# x4 = runif(10)
# y = runif(10)
# testnew = sbm.check_formula(y~x1+x2+x3+x4|x2+x1+x3)
# testnew2 = sbm.check_formula(y~x1+x2|x3)
# zv = testnew$zvars
# xv = testnew$xvars
#Main estimation function for a sbm class of given number of breaks 


#test if we can simplify an argument
