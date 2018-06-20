#' This function computes the cv errors for different smoothing parameter 
#' @param gamma_pool a vector of canadiate smoothing parameter values 
#' @param folds cv folds
#' @param obs observed response, a list 
#' @param timp timepoints corresponding to the observed response, a list 
#' @param basis_intelist a list of integrated B-spline basis function
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param betalist previous DeFEC's coefs
#' @param thres threshold 
#' @param evalc evaluate times 
#' @param maxit maxit 
#' @export select_gamma
#' @import Rsolnp
#' @import fda
#' @import ggplot2
#' @import dplyr
#' @examples
#' \dontrun{
#' }

select_gamma = function(gamma_pool, folds, obs, timep,basis_intelist,spline_basis,thres=1e-4,betalist=NULL,evalc=5,maxit=1e3)
{

cv_res = lapply(gamma_pool , function(gam){
fit_error = c()
for (j in 1:length(folds)){
	test = folds[[j]]
	train = setdiff(1:length(observed),test)
pc_multi= FPCfirst_multi_start(evaltimes=evalc,previous_beta=betalist,obs=obs[train], timep=timep[train],basis_intelist =basis_intelist[train], spline_basis=spline_basis, gamma=gam,threshold=1e-2,minit=20,maxeval=maxit)
pc = pc_multi[[which.min(sapply(pc_multi, function(x) x$options$value))]]
betalist = list()
betalist[[1]]= pc$beta
res = predict_DeFEC(betalist,obs[test], timep[test], basis_intelist =basis_intelist[test],spline_basis,nminus=1)

fit_error  =c(fit_error,res$sigmahat%>%mean)
}
fit_error 
})%>%do.call(rbind,.)

cv_res%>%rowMeans
(gamma_select = gamma_pool[which.min(cv_res%>%rowMeans)])

return(list(gamma_select=gamma_select, cv_res = cv_res, gamma_pool=gamma_pool))

}


FPCfirst_multi_start = function(evaltimes = 3, obs,previous_beta=NULL, timep ,basis_intelist, spline_basis, threshold,minit,gamma, maxeval) {
	 
	results = list()
	 for (j in 1:evaltimes){
     	
	pc_fit = single_FPC (betalist=previous_beta,obs, timep ,basis_intelist,  spline_basis, threshold,minit,gamma,maxeval,randomseed = as.numeric(Sys.time()))
	results[[j]] = pc_fit
	}
	 return(results)
}


