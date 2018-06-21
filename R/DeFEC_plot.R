#' This function plots the estimated Derivatives 
#' @param DeFEC_pred results from function predict_DeFEC
#' @param DeFEC_fd results from function DeFEC
#' @param sample_ids sample id number 
#' @export DeFECplot
#' @import fda
#' @import ggplot2
#' @import dplyr
#' @examples
#' \dontrun{
#' }

DeFECplot = function(DeFEC_pred,DeFEC_fd,sample_ids ) {

observed = DeFEC_pred$observed
timepoints = DeFEC_pred$timepoints

spline_basis = DeFEC_fd$basis
# t = seq(spline_basis$rangeval[1],spline_basis $rangeval[2],len=100)
t = seq(spline_basis$rangeval[1],spline_basis $rangeval[2],len=30)
basis_intemat  = lapply(t, function(x){
	sapply(1:spline_basis$nbasis,function(i){
	basis_integral = function(x){
	fb = function(t) {sapply(t,function(x) {eval.basis(x,spline_basis)[i]%>%as.numeric})}
	integrate(fb,spline_basis$rangeval[1],x)$value
	}
	basis_integral(x)
	})
})%>%do.call(cbind,.)

dpc_integrate  = t(basis_intemat)%*%DeFEC_fd$coefs


plotdata = (dpc_integrate)%*%t(DeFEC_pred$sfit[,-1])
plotdata2 = matrix(DeFEC_pred$sfit[,1]%>%rep(each=length(t)), nrow=length(t))+plotdata
data = reshape2::melt(plotdata2)
names(data) =c('Time',"id","value")
data$Time = t[data$Time]

data_obs = lapply(sample_ids, function(id) {
data.frame(y=observed[[id]], Time=timepoints[[id]], id=id)

})%>%do.call(rbind.data.frame,.)


p1=  ggplot(data%>%filter(id%in%sample_ids),aes(x=Time, y=value,group=id,color=factor(id)))+geom_line(size=1.1,alpha=0.9,col=4)+geom_point(data=data_obs%>%dplyr::filter(id%in%sample_ids),aes(x=Time, y=y, group=id),size=2,col=2)+facet_wrap(~id,scales='free')+theme_bw()+ylab('Y')+xlab('Time')+theme(legend.position = "none",text=element_text(size=rel(4)),legend.text =element_text(size=rel(4)))



dpc=eval.fd(t, DeFEC_fd)


plotdata = (dpc%>%as.matrix)%*%t(DeFEC_pred$sfit[,-1])
data = reshape2::melt(plotdata)
names(data) =c('Time',"id","value")
data$Time = t[data$Time]

p2 = ggplot(data%>%filter(id%in%sample_ids),aes(x=Time, y=value,group=id))+geom_line(size=1.1,col=4)+geom_hline(yintercept=0, linetype=2,col=2)+facet_wrap(~id, scales = "free")+theme_bw()+ylab('Estimated Derivatives')+xlab('Time')+theme(legend.position = "none",text=element_text(size=rel(4)),legend.text =element_text(size=rel(4)))

return(list(fitplot = p1, derivplot=p2))


	}





