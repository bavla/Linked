#
# Pajek.R
# reading and saving of Pajek data in R (matrices and vectors)
#
# by Vladimir Batagelj, November 2019
#
# https://github.com/bavla/cluRC/blob/master/readPajekNet.R
#

clu2vector <- function(f,skip=1){
  read.table(f,skip=skip,colClasses=c("integer"),header=FALSE)$V1
}

vec2vector <- function(f,skip=1){
  read.table(f,skip=skip,colClasses=c("numeric"),header=FALSE)$V1
}

net2matrix <- function(f,skip=0){
# reads a network from Pajek's net file; skip initial comments lines
   L <- readLines(f)
   st <- grep("\\*",L)
   S <- unlist(strsplit(trimws(L[st[1]]),'[[:space:]]+'))
   n <- as.integer(S[2]); n1 <- st[1]+1; n2 <- st[2]-1
   m1 <- st[2]+1; m2 <- length(L); m <- m2-m1+1
   Names <- unlist(strsplit(L[n1:n2],'"'))[2*(1:n)]
   R <- matrix(data=0,nrow=n,ncol=n,dimnames=list(Names,Names))
   S <- unlist(strsplit(trimws(L[m1:m2]),'[[:space:]]+'))
   b <- as.integer(S[3*(1:m)-2]); e <- as.integer(S[3*(1:m)-1]); v <- as.numeric(S[3*(1:m)])
   for(k in 1:m) R[b[k],e[k]] <- R[b[k],e[k]]+v[k]
#   rownames(R) <- Names
   return(R)
}

matrix2net <- function(M,Net="Pajek.net"){
  n <- nrow(M); net <- file(Net,"w")
  cat("% mat2Pajek",date(),"\n*vertices",n,"\n",file=net)
  RN <- row.names(M)
  for(v in 1:n) cat(v,' "',RN[v],'"\n',sep="",file=net)
  cat("*arcs\n",file=net)
  for(v in 1:n) for(u in 1:n) if(M[v,u]!=0) cat(v,u,M[v,u],"\n",file=net)
  close(net)
}

bimatrix2net <- function(M,Net="Pajek.net"){
  n <- nrow(M); m <- ncol(M); net <- file(Net,"w")
  cat("% bip2Pajek",date(),"\n*vertices",n+m,n,"\n",file=net)
  RN <- dimnames(M)[[1]]; CN <- dimnames(M)[[2]];
  for(v in 1:n) cat(v,' "',RN[v],'"\n',sep="",file=net)
  for(u in 1:m) cat(u+n,' "',CN[u],'"\n',sep="",file=net)
  cat("*arcs\n",file=net)
  for(v in 1:n) for(u in 1:m) if(M[v,u]!=0) cat(v,u+n,M[v,u],"\n",file=net)
  close(net)
}

uv2net <- function(u,v,w=NULL,Net="Pajek.net",twomode=FALSE){
  net <- file(Net,"w")
  if(is.null(w)) w <- rep(1,length(u))
  RN <- levels(u); n <- length(RN)
  if(twomode) {CN <- levels(v);  m <- length(CN)}
  U <- as.integer(u); V <- as.integer(v)
  if(twomode) cat("% uv2Pajek",date(),"\n*vertices",n+m,n,"\n",file=net) else
    cat("% uv2Pajek",date(),"\n*vertices",n,"\n",file=net)
  for(i in 1:n) cat(i,' "',RN[i],'"\n',sep="",file=net)
  if(twomode) for(i in 1:m) cat(i+n,' "',CN[i],'"\n',sep="",file=net)
  cat("*arcs\n",file=net)
  for(i in 1:length(u)) cat(U[i],V[i]+twomode*n,w[i],"\n",file=net)
  close(net)
}

vector2clu <- function(C,Clu="Pajek.clu"){
  n <- length(C); clu <- file(Clu,"w")
  cat("% clu2Pajek",date(),"\n*vertices",n,"\n",file=clu)
  cat(C,sep="\n",file=clu)
  close(clu)
}

vector2vec <- function(X,Vec="Pajek.vec"){
  n <- length(X); vec <- file(Vec,"w")
  cat("% vec2Pajek",date(),"\n*vertices",n,"\n",file=vec)
  cat(X,sep="\n",file=vec)
  close(vec)
}
