


normalizar<-function(x){
  return(x/sqrt(sum(x*x)))
}
mds_classic<-function(X,nd=2,plot=TRUE){
  B<-X%*%t(X)
  diagonalizacion<-eigen(B)
  vectores<-as.matrix(apply(diagonalizacion$vectors,2,normalizar))
  valores<-as.numeric(diagonalizacion$values)
  vectores<-vectores[,valores!=0]
  valores<-valores[valores!=0]
  valores<-as.numeric(sort(valores,decreasing = TRUE))
  valores_sel<-as.numeric(valores[1:nd])
  P<-as.matrix(vectores[,order(valores,decreasing = TRUE)[1:nd]],nrow=nrow(X),ncol=nd)
  if(nd>1){
    D<-as.matrix(sqrt(diag(valores_sel)),nd,nd)
    proy<-P%*%D
    if(plot==TRUE&&nd==2){
      plot(proy[,1],proy[,2],asp=1,axes=TRUE)
    }
  }
  else{
    D<-valores_sel[1]
    proy<-P*D
    if(plot==TRUE){
      plot(proy[,1],rep(0,length(proy[,1])),asp=1,axes=TRUE)
    }
  }
  return(proy)
}
cluster_sint<-function(coords){
  #Parámetros generales
  nombre<-coords$nombre
  n<-coords$n
  m<-coords$m
  noisefactorpat_1<-coords$noisefactorpat_1
  noisefactorpat_2<-coords$noisefactorpat_2
  noisefactor<-coords$noisefactor
  margingenes<-coords$margingenes
  marginpacientes<-coords$marginpacientes
  #Parámetros del cluster 1
  mu<-coords$mu
  pc1<-coords$pc1
  ng1<-length(pc1)
  p<-coords$p
  alpha<-coords$alpha
  #Parámetros del cluster 2
  delta<-coords$delta
  pc2<-coords$pc2
  ng2<-length(pc2)
  q<-coords$q
  beta<-coords$beta
  
  ############################################################################
  noisefactorpat_1<-noisefactorpat_1*6
  noisefactorpat_2<-noisefactorpat_2*6
  noisefactor<-noisefactor*6
  
  noisematrix<-matrix(rnorm((m+2*marginpacientes)*(n+2*margingenes),mean=0,sd=noisefactor),m+2*marginpacientes,n+2*margingenes)
  
  
  #h_k:[0,2*pi]---->[0,1] k e {0,...,n1-1}
  
  h_k<-function(s){
    alpha_rand<-rnorm(1,mean=alpha,sd=0.01)
    if((s>=0)&&(s<((1+mu)*pi*(1+alpha_rand)))){
      return(((1+cos(((s/(1+mu))-alpha_rand*pi)))/2)^{1/p})
    } else if((s>=((1+mu)*pi*(1+alpha_rand)))&&(s<=2*pi)){
      return(0)
    }else{
      return("error dominio h_k")
    }
  }
  #g_k:[0,2*pi]---->[0,1] k e {0,...,n2-1}
  g_k<-function(s){
    beta_rand<-rnorm(1,mean=beta,sd=0.01)
    if((s<(1-delta)*pi*(1-beta_rand))&&(s>=0)){
      return(0)
      
    } else if(s>=(1-delta)*pi*(1-beta_rand)&&(s<=2*pi)){
      return(((1+cos(((2*pi-s)/(1+delta))-beta_rand*pi))/2)^(1/q))
      
    }else{
      return("error dominio g_k")
    }
  }
  matrixcdef<-matrix(0,m+2*marginpacientes,n+2*margingenes)+noisematrix
  
  x<-seq(from=0,to=2*pi,length.out = n)
  if(ng1>1){
    cluster1<-function(){
      c1m<-matrix(0,ng1,n)
      for(i in 0:(ng1-1)){
        v<-sapply(x, h_k )
        c1m[i+1,]<-v
      }
      rf<-sample(ng1)
      c1m<-c1m[rf,]
      return(c1m)
    }
    p1<-cluster1()
    #pheatmap(p1,color=colorRampPalette(c("grey","yellow"))(100),cluster_cols = FALSE, cluster_rows = FALSE)
    noisematrixcluster1<-matrix(rnorm(ng1*n,mean=0,sd=noisefactorpat_1),ng1,n)
    p1<-p1+noisematrixcluster1
    #pheatmap(p1,color=colorRampPalette(c("grey","blue"))(100),cluster_cols = FALSE, cluster_rows = FALSE)
    matrixcdef[marginpacientes+pc1,((margingenes+1):(n+margingenes))]<-p1
  }
  
  if(ng2>1){
    cluster2<-function(){
      c2m<-matrix(0,ng2,n)
      for(i in 0:(ng2-1)){
        v<-sapply(x, g_k)
        c2m[i+1,]<-v
      }
      rf<-sample(ng2)
      c2m<-c2m[rf,]
      return(c2m)
    }
    p2<-cluster2()
    #pheatmap(p2,color=colorRampPalette(c("grey","yellow"))(100),cluster_cols = FALSE, cluster_rows = FALSE)
    noisematrixcluster2<-matrix(rnorm(ng2*n,mean=0,sd=noisefactorpat_2),ng2,n)
    p2<-p2+noisematrixcluster2
    #pheatmap(p2,color=colorRampPalette(c("grey","yellow"))(100),cluster_cols = FALSE, cluster_rows = FALSE)
    matrixcdef[(marginpacientes+pc2),((margingenes+1):(margingenes+n))]<-p2
  }
  pci<-intersect(pc1,pc2)
  ngi<-length(pci)
  if(ngi!=0){
    cim<-matrix(0,ngi,n)
    for(i in 0:(ngi-1)){
      random<-runif(1,0.25,0.75)
      v<-(random)*sapply(x, h_k)+(1-random)*sapply(x, g_k)
      cim[i+1,]<-v
    }
    #cim[(0:ngi),(floor(n*0.4):floor(n*0.6))]<-1
    noisematrixclusteri<-matrix(rnorm(ngi*n, mean=0, sd=(noisefactorpat_1+noisefactorpat_2)/2),ngi,n)
    pi<-cim+noisematrixclusteri
    matrixcdef[(marginpacientes+pci),((margingenes+1):(margingenes+n))]<-pi
  }
  matrixf<-as.data.frame(matrixcdef)
  return(matrixf)
}
coords<-list(nombre="X",n=150,m=80,noisefactorpat_1=0.2,noisefactorpat_2=0.2,noisefactor=1,margingenes=0,marginpacientes=0,mu=0,pc1=1:40,p=6,alpha=0.8,delta=0,pc2=41:80,q=6,beta=0.8)
X<-as.matrix(cluster_sint(coords))
proy<-mds_classic(X,nd=2,plot = FALSE)
plot(proy)
plot(proy[coords$pc1,1],proy[coords$pc1,2],col="blue",xlim = c((min(proy[,1])-1),(max(proy[,1])+1)),ylim = c((min(proy[,1])-1),(1+max(proy[,1])))) 
lines(proy[coords$pc2,1],proy[coords$pc2,2],col="red",type = "p")


