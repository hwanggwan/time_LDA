log.hyper <- function(a,b,nii,m,a.a,a.b, b.a,b.b){
  prob <- dgamma(a, a.a, a.b, log = T) + dgamma(b, b.a, b.b,log=T) +
    (m+1)*(log(b)+lgamma(a+b)-lgamma(a))+sum(lgamma(nii+a)-lgamma(nii+1+a+b))
  return(prob)
}


alpha1 <- c(5,1,3,1,1)
alpha2 <- c(1,2,1,4,1)
alpha3 <- c(5,1,1,1,5)

alpha.list <- list(alpha1, alpha2, alpha3)

nt = 30
Nt = nt*3


library(MCMCpack)
phi <- matrix(0, nc=5, nr=Nt)
for(i in 1:3){
  phi[((i-1)*nt+1):(nt*i),] <- rdirichlet(nt, alpha.list[[i]])
}

dn <- rpois(Nt, lambda = 100)
psi <- matrix(0, 5, 25)
for(i in 1:5){
  tmp.alpha <- rep(1,25)
  tmp.alpha[((i-1)*5+1):(i*5)] <- 50
  psi[i,] <- rdirichlet(1, alpha= tmp.alpha)
}


z <- apply(phi, 1, which.max)
#z <- rep(1:5, each=5)

doc.tb <- matrix(0, nc=25, nr=Nt)
colnames(doc.tb) <- letters[1:25]
doc.list <- vector('list', length = Nt)



for(i in 1:Nt){
  zw <- sample(1:5,size=dn[i],replace = T,prob=phi[i,])
  doc <- sapply(1:dn[i], function(v) letters[sample(1:25, size=1, prob=psi[zw[v],])])
  doc.list[[i]] <- doc
  doc.tb[i,] <- table(factor(doc, levels = letters[1:25]))
}

st <- sort(c(1:10,sample(1:10, size=Nt-10, replace=T )))

z.hat <- vector('list', length = Nt)

n.vk <- matrix(0, nr=5, nc=25)
n.mk <- matrix(0, nr=5, nc=Nt)
colnames(n.vk) <- letters[1:25]
for(i in 1:Nt){
  for(w in 1:dn[i]){
    tmp.z <-sample(1:5, size=1)
    z.hat[[i]][w] <- tmp.z
    n.vk[tmp.z,doc.list[[i]][w]] <- n.vk[tmp.z,doc.list[[i]][w]]+1
    n.mk[tmp.z,i] <- n.mk[tmp.z,i] + 1
  }

}
alpha0 <- 1


niter <- 1000
ndisp <- 10
withinIter <- 5

a <- 1
b <- 1

alpha.hat <- t(rdirichlet(Nt,rep(1,5)))
sg <- matrix(0, nr=niter, nc=Nt)
set.seed(123)
for(iter in 1:niter){
  for(iiter in 1:withinIter){
    for(t in 1:Nt){
      for(i in 1:dn[t]){
        w <- doc.list[[t]][i]
        zw <- z.hat[[t]][i]
        n.vk[zw,w] <- n.vk[zw,w]-1
        n.mk[zw,t] <- n.mk[zw,t]-1
        tmp.prob <- log((n.vk[,w]+1))-log(apply(n.vk,1,sum)+25) +log(n.mk[,t]+alpha.hat[,st[t]])
        tmp.z <-sample(1:5, size=1, prob=exp(tmp.prob - max(tmp.prob)))
        z.hat[[t]][i] <- tmp.z
        
        n.vk[tmp.z,w] <- n.vk[tmp.z,w] + 1
        n.mk[tmp.z,t] <- n.mk[tmp.z,t] + 1
        #      n.vk <- sapply(letters[1:25], function(w) table(factor(z.hat[unlist(doc.list)==w], levels = 1:5)))
        #      n.mk <- sapply(1:Nt, function(m) table(factor(z.hat[(cum.dn[m]+1):cum.dn[m+1]],levels = 1:5)))
      }
    }
    
  }

  st.set <- unique(st)
  phi.hat <- sapply(1:Nt, function(t){
    tmp.param <- n.mk[,t]+alpha.hat[st[t]]
    tmp.phi <- rdirichlet(1, tmp.param) + 1e-08
    tmp.phi <- tmp.phi/sum(tmp.phi)
    return(tmp.phi)})
#  phi.hat <- apply(n.mk, 2, function(x) x/sum(x))
  alpha.hat <- sapply(st.set, function(s){
    est.dir <- dirichlet::fit.dirichlet(t(phi.hat[,st==s]))
    tmp.k <- est.dir$most.likely.k
    return(est.dir$p*ifelse(is.na(tmp.k),1,tmp.k))
    }  )
  #alpha.hat <- sapply(st.set, function(s) rowSums(n.mk[,st==s, drop=F]) + alpha0)
#  alpha.mode <- apply(alpha.hat, 2, function(x) rdirichlet(1, x))
  alpha.mode <- alpha.hat
  
  
    if(st[1]!=st[2]){
      prob <- c(log(a*b/(b+a)^2)+log(ddirichlet(phi.hat[,1], alpha=alpha.mode[,st[1]])),
                log(b)+log(a+(sum(st[2:(Nt-1)]==st[2])-1)) -log(a+b) -
                  log(a+(sum(st[2:(Nt-1)]==st[2])-1)+b) + log(ddirichlet(phi.hat[,1], alpha=alpha.mode[,st[2]])))
      prob <- exp(prob-max(prob))
      st[1] <- sample(st[1:2], size=1, prob=prob)
    }
    
    
    for(t in 2:(Nt-1)){
      tmp.st <- st[c(t-1,t+1)]
      if(diff(tmp.st)==0) next
      n11 <- sum(st[1:(t-1)]==tmp.st[1])-1
      n22 <- sum(st[(t+1):(Nt-1)]==tmp.st[2])-1
      p1 <- log(n11+a)+log(b) -log(n11+b+a)-log(n11+1+b+a) + log(ddirichlet(phi.hat[,t], alpha=alpha.mode[,tmp.st[1]]))
      p2 <- log(b) + log(n22+a) - log(n11+b+a) - log(n22+b+a) + log(ddirichlet(phi.hat[,t], alpha=alpha.mode[,tmp.st[2]]))
      prob = exp(c(p1,p2)-max(c(p1,p2)))
      st[t] <- sample(tmp.st, size=1, prob=prob)
    }
    
    if(st[Nt]!=st[Nt-1]){
      n33 <- sum(st[1:(Nt-1)]==st[Nt-1])-1
      p1 = log(n33+a)-log(n33+b+a)+log(ddirichlet(phi.hat[,Nt], alpha=alpha.mode[,st[Nt-1]]))
      p2 = log(b)-log(n33+b+a) + log(ddirichlet(phi.hat[,Nt], alpha=alpha.mode[,st[Nt]]))
      prob=exp(c(p1,p2)-max(c(p1,p2)))
      st[Nt] <- sample(c(st[Nt-1],st[Nt]), size=1, prob=prob)
    }
  
  st <- as.integer(factor(st))
  sg[iter,] <- st
  st.set <- unique(st)
  nii <- sapply(st.set, function(k) sum(st==k)-1)
  m <- length(st.set)-1
  
  a.new <- truncnorm::rtruncnorm(1,a=0, mean=a)
  testp <- log.hyper(a.new,b,nii,m,1,1,1,1) - log.hyper(a,b,nii,m,1,1,1,1) + log(pnorm(a)) - log(pnorm(a.new))
  b.new <- truncnorm::rtruncnorm(1,a=0, mean=b)
  testp <- log.hyper(a,b.new,nii,m,1,1,1,1) - log.hyper(a,b,nii,m,1,1,1,1) + log(pnorm(b)) - log(pnorm(b.new))
  if(testp < log(runif(1))) b <- b.new
  
  
  if(iter %% ndisp ==0){
    print(iter)
    print('phi.hat')
    print(apply(n.mk, 2, which.max))
    print('st')
    print(st)
    }
}

cp <- sg[160:380,]
plot(rowMeans(apply(cp, 1, function(x) c(F,diff(x)!=0))), type='h')
abline(v=c(31,61), col='blue', lty=2, lwd=2)

alpha.hat

aricode::ARI(z,apply(phi.hat, 2, which.max))
aricode::NMI(z,apply(phi.hat, 2, which.max))

phi.hat
