aiFun<-function(model, Dimnames=NULL){
  AI.list<-model$ai
  k<-length(model$gammas)
  tmp.index<-matrix(c(1,2,4,7,11,16,22,29,37,46,
	              2,3,5,8,12,17,23,30,38,47,
  	              4,5,6,9,13,18,24,31,39,48,
	              7,8,9,10,14,19,25,32,40,49,
  	              11,12,13,14,15,20,26,33,41,50,
	              16,17,18,19,20,21,27,34,42,51,
    	              22,23,24,25,25,27,28,35,43,52,
	              29,30,31,32,33,34,35,36,44,53,
		      37,38,39,40,41,42,43,44,45,54,
		      46,47,48,49,50,51,52,53,54,55), nrow=10, ncol=10)

  index<-as.vector(tmp.index[1:k, 1:k])

  AI.cov<-matrix(rep(AI.list,2)[index], nrow=k, ncol=k)
  AI.cor<-cov2cor(AI.cov)
  AI.mat<-matrix(0, nrow=k, ncol=k)
  AI.mat[upper.tri(AI.cor, diag=FALSE)]<-AI.cor[upper.tri(AI.cor, diag=FALSE)]
  AI.mat[lower.tri(AI.cov, diag=TRUE)]<-AI.cov[lower.tri(AI.cov, diag=TRUE)]
  if(!is.null(Dimnames)) dimnames(AI.mat)<-list(Dimnames, Dimnames)
return(round(AI.mat, 5))
}

