#' @title P values calculation in The multivariable Cox model Using LRT
#'
#' @description The function creates a p-value table for multivariable Cox Proportional Hazard regression model.
#' @param time A variable for time to event
#' @param event A variable for event indicator (1 indicates event, 0 is censored)
#' @param covariate A vector of variables included in the multivariable analysis
#' @param data A data set
#' @export
#' @import survival
#' @return A matrix with variables and p values in the multivariable Cox model
#' @examples \dontrun{
#' library(survival)
#' df <- veteran
#' res <- lr.surv.pval(time="time", event="status",
#'             covariate=c("age","prior","trt","celltype","karno","diagtime"),
#'             data=df)
#' }


lr.surv.pval <- function(time, event, covariate, data){
  # Start function
  df <- data
  cov.num   <- length(covariate)
  select    <- c(time,event,covariate)

  tab <- NULL
  for (i in 1:cov.num){
    temp <- df[,select]
    temp <- temp[complete.cases(temp), ]

    cov.f <- paste(covariate, collapse="+")
    cov.r <- paste(covariate[-i], collapse="+")
    out   <- paste("Surv(",time,",",event,")", "~",sep='')

    mod.f <- as.formula(paste(out, cov.f))
    mod.r <- as.formula(paste(out, cov.r))

    model.f <- coxph(mod.f, data=temp)
    model.r <- coxph(mod.r, data=temp)

    tmp <- anova(model.f, model.r)$`P(>|Chi|)`[2]
    pval<- round(tmp,10)
    res <- c(covariate[i],pval)
    tab <- rbind(tab,res)
    rownames(tab) <- NULL
    colnames(tab) <- c("variable","pval")
  }
  tab <- as.data.frame(tab)
  tab[,2] <- as.numeric(as.character(tab[,2]))
  return(tab)
}






#' @title A backward elimination model building in Cox regression
#'
#' @description The function builds or finalizes a multivariable Cox regression model based on a LRT p-value.
#' @param time A variable for time to event
#' @param event A variable for event indicator (1 indicates event, 0 is censored)
#' @param covariate A vector of variables included in the multivariable analysis
#' @param sig A significance level that allows the variable in the finalized model
#' @param data A data set
#' @export
#' @import survival
#' @return A matrix with variables and p values in the multivariable Cox model
#' @examples \dontrun{
#' library(survival)
#' df <- veteran
#' res <- surv.pval(time="time", event="status",
#'             covariate=c("age","prior","trt","celltype","karno","diagtime"),
#'             sig=0.05, data=df)
#' }

surv.pval <- function(time, event, covariate, sig=0.05, data){

  len <- length(covariate)
  res <- mod.building::lr.surv.pval(time, event, covariate, data)
  significance <- sig

  while (len>2 & len != sum(res$pval<significance)) {

    covariate <- covariate[-which.max(res$pval)]


    len <- length(covariate)
    res <- mod.building::lr.surv.pval(time, event, covariate, data);
  }

  if (len != sum(res$pval<significance)) {
    print("There is no mutivariable fit.")
    print(res)
  } else {
    return(res)
  }
}









