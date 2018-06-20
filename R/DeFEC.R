#' This function computes the first K DeFEC
#' @param observed a list of observed longitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param timepoints a list of observed time points 
#' @param K the number of DeFEC
#' @param basis_intematlist a list of integrated B-spline basis function
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param gammas a vector of length K, the smoothing parameter 
#' @export DeFEC
#' @import fda
#' @import dplyr
#' @import Rsolnp
#' @examples
#' \dontrun{
#' 
#' }

DeFEC = function(observed, timepoints, K=3,spline_basis, basis_intematlist ,thresh=1e-3, gammas=NULL, maxit=5e2)
	{

	if (is.null(gammas)) gammas = rep(0,K)
	for (p in 1:K){

	if (p==1) betalistp = NULL

	pc_fit = single_FPC(betalist= betalistp, obs=observed,timep=timepoints,basis_intelist=basis_intematlist, spline_basis=spline_basis,threshold=thresh,minit=20,gamma=gammas[p],maxeval=maxit)

	betalistp[[p]]  = pc_fit$beta
	tempbetalist = NULL
	tempbetalist[[1]] = pc_fit$beta
	fit_pc = pred_FOB(betalist=tempbetalist,ylist=observed, tlist=timepoints,basis_intelist=basis_intematlist,spline_basis=spline_basis, nminus=1,intercept=TRUE)
	observed = fit_pc$residuals

	}
	return(betalistp)
	

}