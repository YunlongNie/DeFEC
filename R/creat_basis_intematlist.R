#' This function create the integrated bspline functions
#' @param timepoints a list of observed time points
#' @param spline_basis the B-spline basis for the fitted FEC
#' @export create_basis_intematlist
#' @import fda
#' @import dplyr
#' @import Rsolnp
#' @examples
#' \dontrun{
#' 
#' }

create_basis_intematlist = function(timepoints, spline_basis){
	basis_intematlist = lapply(timepoints,function(t){
	basis_intemat  = lapply(t, function(x){
		sapply(1:spline_basis$nbasis,function(i){
		basis_integral = function(x){
		fb = function(t) {sapply(t,function(x) {eval.basis(x,spline_basis)[i]%>%as.numeric})}
		integrate(fb,spline_basis$rangeval[1],x)$value
		}
		basis_integral(x)
		})
	})%>%do.call(cbind,.)
	basis_intemat  
	})
	return(basis_intematlist)
	}

