source("misc.R")
library(survival)
#set.seed(134)

 B <-200
 N <- 100
 hatb.l <- matrix(NA, nrow=B, ncol=2)
 hatb <- rep(NA, B) 
 mm <- yy<- rep(NA, B) 
for( bb in 1: B )
{
  print(bb)
  
  
  mydata <- gendata(N= N , a=8, b=1, r=0, k=1, tau=10 ) 
  redata <- data.frame( mydata$re )  
  lwdata <- data.frame( mydata$lw )  
    
  #--------  Estimation using LWYY method ----------#

   fit.l <- summary( coxph( Surv(start, stop, event)~x+cluster(id), 
                           data= lwdata, robust=TRUE))
   hatb.l[bb,] <- fit.l$coef[1,c(1,4)]
   
 #--------  Estimation with pairwise pseudolikelihood ----------#
 
  idfreq <- table(redata$id)
  uid <- unique(redata$id)
  mi <- tapply( redata$t, redata$id, function(x) length(x[x!=0]) )
  xi <- tapply( redata$x, redata$id, mean)
  yi <- tapply( redata$y, redata$id, mean)
  #wi <- tapply( redata$w, redata$id, mean)
  mm[bb]<- mean(mi)
  yy[bb]<- mean(yi)

  lastrow  <- cumsum( idfreq )
  firstrow <- c(1, 1+lastrow[-N] )

 #--------  calculate R ----------#
  
  R <- NULL
  ct <-0
  for( i in 1: (N-1) ) {
    if( firstrow[i]!=lastrow[i] | redata$t[firstrow[i]]!=0 ){
      datai <- redata[ firstrow[i]:lastrow[i], ]      

      for ( k in (i+1): N ){ 
        if ( firstrow[k]!= lastrow[k]| redata$t[firstrow[k]]!=0 )
        { 
          ct <- ct+1  
          datak <- redata[ firstrow[k]:lastrow[k], ]    
          R <- rbind( R, cbind(count=ct, i=i,k=k,r=c( R.cal(datai, datak) ) ))
        }     
        datak <- yik <- ti <-tk <- NULL
      }
    }
    datai <- NULL
  }

  R <- data.frame(R)
  R <- R[R$r!=-99999,]

  lik <- function(b, R=R$r) {
           sum( log( 1+ exp(-R*b)) )    }

  hatb[bb] <- nlm( lik, p=1, R=R$r )$estimate

 #--------  Calculate asymptotic variance ----------#

   ### Asymptotic variance ###
   fith <- data.frame(   hik(hatb[bb], R)  ) 
   tmph <- data.frame( 
             cbind( i=rep(1:N, rep(N,N)), k=rep(1:N, N), h=0, dh=0) )
   tmph[ (fith$i-1)*N+fith$k, 3:4] <- fith[,3:4] 
   tmph[ (fith$k-1)*N+fith$i, 3:4] <- fith[,3:4] 
   tmph <- tmph[ tmph$i!=tmph$k, ]

   sig <- 4 * sum( tapply( tmph$h, tmph$i,  function(x){  
                       sum(outer(x, x, "*")) }    ))/N/(N-1)^2 
           #- 
           #4 * sum( fith$h )/choose(N,2)

    v2 <- sum(fith$dh)/choose(N,2)
 
   asymvar[bb] <- sig/v2^2/N
 

   ### Wald test ###

   fith0 <- data.frame(  hik( r0, R))    
   tmph <- data.frame( 
             cbind( i=rep(1:N, rep(N,N)), k=rep(1:N, N), h=0, dh=0) )
   tmph[ (fith0$i-1)*N+fith0$k, 3:4] <- fith0[,3:4] 
   tmph[ (fith0$k-1)*N+fith0$i, 3:4] <- fith0[,3:4] 
   tmph <- tmph[ tmph$i!=tmph$k, ]


   sig0 <- 4 * sum( tapply( tmph$h, tmph$i,  function(x){  
                       sum(outer(x, x, "*")) }    ))/N/(N-1)^2 - 
           4 * sum( fith0$h )/choose(N,2)
   v20 <- sum(fith0$dh)/choose(N,2)
 
 
   #wald1[bb]<- N*(hatb[bb]-r0)^2*v2^2/sig0
   wald2[bb]<- N*(hatb[bb]-r0)^2*v20^2/sig0


 }
 
   mean(mm)
   mean(hatb.l[,1])
   sd(hatb.l[,1])
   sqrt(mean(hatb.l[,2]^2) )
   mean(hatb )



   #----------  estimate lambda_0 ------------#
   
   #----------   estimate gamma   -------------#

   #----------fit LWYY model -----------#
  lwyy <- data.frame( sim$lwyy )
  l.fit <- summary( coxph( Surv(start, stop, event)~Xt+X, data=lwyy ))  
  lhatb[bb,] <- l.fit$coef[1,c(1,3,5)]
  lhatr[bb,] <- l.fit$coef[2,c(1,3,5)]
  lcib[bb,] <- log( l.fit$conf.int[1,3:4])
  lcir[bb,] <- log( l.fit$conf.int[2,3:4])
}

  cat(file="simout.txt",
      "\n-------------------------------------------------\n",
      "N =", N, ", r =", beta, ", b =", gamma, 
      ", m =", round(mean(mm),1), ", y =", round(mean(yy),1),"\n",
      "Proposed:\n",
      "\t bias (r)   =", round( mean(hatb)-beta, 3), "\n", 
      "\t emp sd(r)  =", round( sd(hatb), 3), "\n",
      "\t asym sd(r) =", round( mean(asymvar), 3), "\n",
      "\t cov pr(r)  =", round( mean(beta>hatb-1.96*sqrt(asymvar) 
                               & beta< hatb+1.96*sqrt(asymvar)), 2), "\n",
      "\t rejec pr   =", round( mean(pchisq(wald2,1)>0.95), 2), "\n",
      "LWYY:\n",
      "\t bias (r)   =", round( mean(lhatb[,1])-beta, 3), "\n",
      "\t emp sd(r)  =", round( sd(lhatb[,1]), 3), "\n",
      "\t asym sd(r) =", round( mean(lhatb[,2]), 3), "\n",
      "\t cov pr(r)  =", round( mean( lcib[,1]<beta & beta<lcib[,2]), 2), "\n",
      "\t rejec pr   =", round( mean(lhatb[,3]<0.05), 2), "\n\n",

      "\t bias (b)   =", round( mean(lhatr[,1])-gamma, 3), "\n",
      "\t emp sd(b)  =", round( sd(lhatr[,1]), 3), "\n",
      "\t asym sd(b) =", round( mean(lhatr[,2]), 3), "\n",
      "\t cov pr(b)  =", round( mean( lcir[,1]<beta & beta<lcir[,2]), 2), "\n",
      "\t rejec pr   =", round( mean(lhatr[,3]<0.05), 2), "\n",
      "---------------------------------------------------\n")
      
    
  
 