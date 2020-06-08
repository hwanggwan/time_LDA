log.hyper <- function(a,b,nii,m,a.a,a.b, b.a,b.b){
  prob <- dgamma(a, a.a, a.b, log = T) + dgamma(b, b.a, b.b,log=T) +
    (m+1)*(log(b)+lgamma(a+b)-lgamma(a))+sum(lgamma(nii+a)-lgamma(nii+1+a+b))
  return(prob)
}

sample_conparam <- function(curr_param,numdata,tot_class,aa,bb,numiter){
  
  n_group<-length(numdata)
  alpha <- curr_param
  for(i in 1:numiter){
    wvec_sum = 0
    svec_sum = 0
    
    for(j in 1:n_group){
      if(numdata[j]>0){
        wvec_sum <- wvec_sum + log(rbeta(1,alpha+1,numdata[j]))  ### beta
        if(runif(1)<numdata[j]/(alpha+numdata[j])) svec_sum <- svec_sum + 1 ###bernoulli
      }
      
    }
    gammaa <- aa + tot_class - svec_sum
    gammab <- bb - wvec_sum
    alpha = rgamma(1,gammaa,gammab)
  }
  
  return(alpha)
}

'update.s' <- function(yobs, py,  P){
  
  ns <- nrow(P)
  n <- ncol(yobs)
  py <- t(py)  
  trans  <- matrix(0,ns,ns)
  nstate <- matrix(0,ns,1)
  M <- matrix(0,ns,n)
  pr1 <- matrix(1,ns,1)
  pstyt1 <- matrix(0,ns, 1)
  logP <- log(P)
  
  for( t in n:1){
    yt <- yobs[,t]
    
    
    unnorm_pstyt <- numeric(ns)
    pstyt <- numeric(ns)
    
    if(t == n){
      pstyt1 <- py[,t]
    }else{
      pstyt1 <- M[,(t+1)]+py[,t] 
    }
    
    for(j in 1:ns){
      Phold <- t(logP[j,])+pstyt1
      unnorm_pstyt[j] <- log(sum(exp(Phold-max(Phold))))+max(Phold)
    }
    M[,t]<-unnorm_pstyt
  }
  
  st = 1
  s <- numeric(n)
  ps <- matrix(0,n,ns)
  pstyn <- numeric(ns)
  Pst_1 <- numeric(ns)
  unnorm_pstyn <- numeric(ns)
  cump <-numeric(ns)
  durs <-numeric(n)
  s_norep <-numeric(n)
  
  for(t in 1:n){
    if(t==1){
      unnorm_pstyn <- M[,t+1]+py[,t]
    }else{
      st = s[t-1]
      Pst_1 <- t(logP[st,])
    }
    
    if(t == n){
      unnorm_pstyn <- py[,t] + Pst_1
    }else{
      unnorm_pstyn <- M[,t+1] + py[,t] + Pst_1
    }
    pstyn <- unnorm_pstyn - log(sum(exp(unnorm_pstyn - max(unnorm_pstyn))))- max(unnorm_pstyn)
    pstyn <- exp(pstyn)
    s[t] <- sample(1:ns,size=1,prob=pstyn)
    
    if(t > 1){
      trans[st,s[t]] <- trans[st,s[t]]+1
    }
    
    nstate[s[t]] <- nstate[s[t]]+1
    ps[t,] <- pstyn      
    
  }
  
  result<-list(s=s,ps=ps,trans=trans,nstate=nstate)
  
  return(result)
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

#st <- sort(c(1:10,sample(1:10, size=Nt-10, replace=T )))


alpha0 <- 1


#niter <- 1000
niter <- 1
ndisp <- 10
withinIter <- 5

a <- 1
b <- 1

a.alpha = 1
b.alpha = 0.01
a.gamma = 1
b.gamma = 0.01
a.theta = 100
b.theta = 1
K = 10

st <- sample(1:10, size=Nt, replace = T)

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
alpha_p_kappa <- rgamma(1,a.alpha,b.alpha)
sticky.theta <- rbeta(1,a.theta,b.theta)
sticky.alpha <- (1-sticky.theta)*alpha_p_kappa
sticky.kappa <- sticky.theta * alpha_p_kappa
sticky.gamma <- rgamma(1,a.gamma,b.gamma)
sticky.gamma_prime <- numeric(K)


#alpha.hat <- t(rdirichlet(Nt,rep(1,5)))
alpha.hat <- t(rdirichlet(K,rep(1,5)))



P <- matrix((1 - sticky.theta) / (K - 1), nrow = K, ncol = K)
diag(P) <- sticky.theta


niter <- 200
sg <- matrix(0, nr=niter, nc=Nt)
set.seed(123)

for(iter in 1:niter){
  
  #Sampling for word topic
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
      }
    }
    
  }
  

  #Estimating alpha
  
  phi.hat <- sapply(1:Nt, function(t){
    tmp.param <- n.mk[,t]+alpha.hat[st[t]]
    tmp.phi <- rdirichlet(1, tmp.param) + 1e-08
    tmp.phi <- tmp.phi/sum(tmp.phi)
    return(tmp.phi)})

  st.set <- unique(st)
  for(s in st.set){
    est.dir <- dirichlet::fit.dirichlet(t(phi.hat[,st==s]))
    tmp.k <- est.dir$most.likely.k
    alpha.hat[,s] <- est.dir$p*ifelse(is.na(tmp.k),1,tmp.k)
  }
  ind.non <- (1:K)[!(1:K %in% st.set)]
  alpha.hat[,ind.non] <- t(rdirichlet(length(ind.non), rep(1,5)))

  alpha.mode <- alpha.hat
  
  
  #Sampling St
  
  py <- apply(alpha.mode, 2, function(x) log(ddirichlet(t(phi.hat),x)))
  s.out <- update.s(phi.hat, py, P)
  st <- s.out$s
  sg[iter,] <- st
  
  #### update for dishes
  
  rest_dishes <- matrix(0,K,K)
  rest_dishes_over <- matrix(0, K, K) 
  sum_w <- matrix(0, K, 1)
  
  for(j in 1:K){
    for(k in 1:K){
      n_jk <- s.out$trans[j,k]
      if(n_jk==0){
        rest_dishes[j,k] <- 0
      }else{
        m_jk = 1
        for(h in 1:n_jk){
          m_num <- sticky.alpha*sticky.gamma_prime[k]
          if(j==k) m_num <- m_num + sticky.kappa
          m_jk <- m_jk + (runif(1) < m_num/(m_num+h))
          
        }
        rest_dishes[j,k] <- m_jk
        rest_dishes_over[j,k] <- m_jk
      }
    }
    if(rest_dishes[j,j]>0){
      gam_thet = sticky.gamma_prime[j] * (1-sticky.theta)
      pp <- sticky.theta/(gam_thet + sticky.theta)
      sum_w[j] <- rbinom(1,rest_dishes[j,j],pp)
      rest_dishes_over[j,j] <- rest_dishes[j,j] - sum_w[j]
    }
    
  }
  
  #end of updating dish counts
  
  #update gamma_prime
  gamma_prime_dir <- sticky.gamma/K + t(colSums(rest_dishes_over))
  gamma_prime <- rdirichlet(1,gamma_prime_dir)
  
  
  #update P
  
  for(j in 1:K){
    p_dirich_params<-matrix(K,1)
    for(i in 1:K){
      p_dirich_params[i] <- sticky.alpha*gamma_prime[i]+s.out$trans[j,i]
      if(i == j) p_dirich_params[i] <- p_dirich_params[i] + sticky.kappa
    }
    
    if( min(p_dirich_params+1) <= 1) p_dirich_params <- p_dirich_params + .Machine$double.eps
    
    P[j,] <- rdirichlet(1,p_dirich_params)
  }
  
  #update DP concentration parameters
  
  Nkdot <- rowSums(s.out$trans)
  Mkdot <- rowSums(rest_dishes)
  Nkdot_valid <- Nkdot[Nkdot>0]
  Mkdot_valid <- Mkdot[Mkdot>0]
  
  ak0  = alpha_p_kappa
  gamma0 = sticky.gamma
  Mdotk <- (colSums(rest_dishes_over)>0)
  Kbar <- sum(Mdotk)
  Mbar_tot <- sum(rest_dishes_over)
  
  #update alpha+kappa
  alpha_p_kappa <- sample_conparam(curr_param = ak0, numdata = Nkdot_valid,
                                   tot_class = sum(Mkdot_valid),
                                   aa = a.alpha,
                                   bb = b.alpha,
                                   numiter = 50)
  #update sticky.gamma
  sticky.gamma <- sample_conparam(curr_param = gamma0,
                                  numdata = Mbar_tot,
                                  tot_class = Kbar,
                                  aa = a.gamma,
                                  bb = b.gamma,
                                  numiter = 50)
  #update sticky.theta
  sticky.theta <- rbeta(1, a.theta + sum(sum_w), b.theta + (sum(rest_dishes)-sum(sum_w)) )
  
  sticky.kappa <- sticky.theta*alpha_p_kappa
  sticky.alpha <- (1-sticky.theta) * alpha_p_kappa
  
  
  if(iter %% ndisp ==0){
    print(iter)
    print('phi.hat')
    print(apply(n.mk, 2, which.max))
    print('st')
    print(st)
  }
}

cp <- sg
plot(rowMeans(apply(cp, 1, function(x) c(F,diff(x)!=0))), type='h')
abline(v=c(31,61), col='blue', lty=2, lwd=2)

aricode::ARI(z,apply(phi.hat, 2, which.max))
aricode::NMI(z,apply(phi.hat, 2, which.max))
