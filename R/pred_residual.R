#' This function compute the residuals
#' @param ylist a list of observed longitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param tlist a list of observed time points
#' @param betalist the coeffcients of the DeFEC
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param basis_intelist a list of integrated B-spline basis function
#' @import fda
#' @import dplyr
#' @examples
#' \dontrun{
#' 
#' }


pred_FOB = function(betalist,ylist, tlist,basis_intelist,spline_basis, nminus=1,intercept=TRUE){

cfits = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
	 colSums(rep(betalist[[i]],each=length(timei))%>%matrix(nrow=spline_basis$nbasis,byrow = T)*basis_intelist[[x]])
	})%>%do.call(cbind,.)

	index_pc = min(length(ylist[[x]])-nminus, length(betalist))
	if(index_pc<=0) 
		{return(rep(0,length(betalist)+1))}
	 else {
	
	# lm(ylist[[x]]~xmat[,1:1])%>%coef%>%as.numeric
	if(intercept) {
		cfit = lm(ylist[[x]]~xmat[,1:index_pc])%>%coef%>%as.numeric
	} else {
		cfit = lm(ylist[[x]]~0+xmat[,index_pc])%>%coef%>%as.numeric
	}


	# if(length(cfit)<length(betalist)+1) cfit = c(cfit, rep(0,length(betalist)+1-length(cfit)))
	return(cfit)
	}	
})%>%do.call(rbind,.)

# yfits_de= lapply(1:nrow(cfits), function(i){
# 	cfit = cfits[i,-1]
# 	rowSums(mapply("*",cfit, betalist,SIMPLIFY=TRUE))
# }
# )%>%do.call(cbind,.)

# yfitsfd_de = fd(yfits_de, spline_basis)


residuals = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
	 colSums(rep(betalist[[i]],each=length(timei))%>%matrix(nrow=spline_basis$nbasis,byrow = T)*basis_intelist[[x]])
	})%>%do.call(cbind,.)

	index_pc = min(length(ylist[[x]])-nminus, length(betalist))
	if(index_pc<=0) 
		{return(rep(0,length(ylist)))}
	 else {
		if(intercept) {
			cfit = lm(ylist[[x]]~xmat[,1:index_pc])%>%residuals%>%as.numeric
		} else {
			cfit = lm(ylist[[x]]~0+xmat[,index_pc])%>%residuals%>%as.numeric
		}

	return(cfit)
	}	
})

sigmahat = sapply(residuals, function(x) mean(x^2))
list(sigmahat = sigmahat,residuals=residuals,sfit = cfits)

}
