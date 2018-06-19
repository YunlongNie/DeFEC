#' This function computes the first DeFEC
#' @param obs a list of observed longitudinal data; each element in the list is a numeric vector for one subject's observed data
#' @param timep a list of observed time points
#' @param betalist a list of previous DeFEC's coefficients, the resulting one should be orthogonal to those previous computed ones. 
#' @param basis_intelist a list of integrated B-spline basis function
#' @param spline_basis the B-spline basis for the fitted FEC
#' @param gamma a positive number, the smoothing parameter 
#' @export single_FPC
#' @import fda
#' @import dplyr
#' @import Rsolnp
#' @examples
#' \dontrun{
#' 
#' }



single_FPC  = function(betalist=NULL, obs,timep, basis_intelist, spline_basis,threshold=1e-3,minit=20,gamma=0,maxeval=1e3){
	beta1 = rnorm(spline_basis$nbasis)
	library(Rsolnp)
	thresh = 1
	it = 1
	beta1_before = rep(0, length(beta1))
	value = -1e14
	pc_fit = fd(beta1, spline_basis)
	pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
	beta1 = coef(pc_fit)%>%as.numeric
	R = inprod(spline_basis, spline_basis,2,2)
	E = inprod(spline_basis,spline_basis,0,0)
	beta_all = list()
	value_all = c()
	while (thresh >  threshold|it<minit){
	beta1_before = beta1
	value_before = value

	alpha_fit = function(i){
		print(i)
		xmat = colSums(rep(beta1,each=length(timep[[i]]))%>%matrix(nrow=spline_basis$nbasis,byrow = T)*basis_intelist[[i]])


		if(is.null(betalist)) {
			lm_fit = lm(obs[[i]]~ xmat)
		} else {
			lm_fit = lm(obs[[i]]~ 0+xmat)
		}
		
		lm_fit%>%coef%>%as.numeric%>%return
	}

	sfit= lapply(1:length(obs),alpha_fit)%>%do.call(rbind,.)

	residual_fit = function(i){
		xmat = colSums(rep(beta1,each=length(timep[[i]]))%>%matrix(nrow=spline_basis$nbasis,byrow = T)*basis_intelist[[i]])
		if(is.null(betalist)) {
			lm_fit = lm(obs[[i]]~ xmat)
		} else {
			lm_fit = lm(obs[[i]]~ 0+xmat)
		}
		
		(lm_fit%>%residuals)^2

	
	}

	rfit= lapply(1:length(obs),residual_fit)%>%do.call(c,.)
	value = mean(rfit)+gamma*inprod(pc_fit,pc_fit,2,2)/inprod(pc_fit,pc_fit)
	value_all = c(value_all,value)


	yem = lapply(1:length(obs),function(i){
		obs[[i]] - sfit[i,1]
		})%>%do.call(c,.)

	sfit_last = sfit[,ncol(sfit)]

	xalpha = lapply(1:length(obs),function(i){

		t(sfit_last[i]*basis_intelist[[i]])

	})%>%do.call(rbind,.)

	if (gamma==0&is.null(betalist))
	{

	beta1 = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)), t(xalpha%>%as.matrix)%*%as.matrix(yem))
	beta1 = beta1/sqrt(as.numeric(t(beta1)%*%E%*%beta1))
	loss = mean((yem  - xalpha%*%beta1)^2)

	} else {

		fn1=function(x)
		{
		mean((yem  - xalpha%*%x)^2)+gamma*as.numeric(t(x)%*%R%*%x)/as.numeric(t(x)%*%E%*%x)
		}
		if(length(betalist) >0)
			{
			betamat=  do.call(rbind,betalist)

			eqn1=function(x){
			z = as.numeric(betamat%*%E%*%x)
			return(z)
			}

			value_now = fn1(beta1)
			loss = value_now - gamma*as.numeric(t(beta1)%*%R%*%beta1)/as.numeric(t(beta1)%*%E%*%beta1)
		
			res_aug = solnp(beta1, fun = fn1, eqfun = eqn1, eqB = rep(0,length(betalist)),control=list(trace=0))


			} else {
				value_now = fn1(beta1)
				loss = value_now - gamma*as.numeric(t(beta1)%*%R%*%beta1)/as.numeric(t(beta1)%*%E%*%beta1)
				
				res_aug = solnp(beta1, fun = fn1,control=list(trace=0))
				}

		it_inside = 1
		b1 = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)), t(xalpha%>%as.matrix)%*%as.matrix(yem))
		b1 = b1/sqrt(as.numeric(t(b1)%*%E%*%b1))
		while(res_aug$values[length(res_aug$values)] >= value&it>2)
			{

		 
			if (length(betalist)>0){
			
				res_aug = solnp(b1+rnorm(length(b1)), fun = fn1, eqfun = eqn1, eqB = rep(0,length(betalist)),control=list(trace=0))

			} else {
				
				res_aug = solnp(b1+rnorm(length(b1)), fun = fn1,control=list(trace=0))
			}
			
			it_inside = it_inside + 1
			if(it_inside>20) break
	
			}

		if(res_aug$values[length(res_aug$values)] < value) {

			beta1 = res_aug$par
		}

		}

	beta_all[[it]] = beta1
	pc_fit = fd(beta1, spline_basis)
	if(inprod(pc_fit)<0) {beta1 = -beta1;pc_fit = fd(beta1, spline_basis)}
	beta1 = coef(pc_fit)%>%as.numeric
	thresh =  max(abs(beta1_before-beta1))
	it = it+1
	if(it%%100==0) {print(it);print(as.numeric(thresh));print(as.numeric(value))}
	if(it>maxeval) break
	}
	pc_fit = fd(beta1, spline_basis)
	pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
	beta1 = coef(pc_fit)%>%as.numeric


	return(list(beta=beta1,pc_fit = pc_fit,sfit=sfit, thresh = thresh,it=it,value=value,gamma=gamma,beta_all = beta_all,value_all = value_all,previous_beta= betalist))
}

