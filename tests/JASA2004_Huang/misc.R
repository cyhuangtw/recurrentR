



#------   Setup for data generation  ------#  
#   lambda_0(t) = exp(kt)                  # 
#   phi(t; x) = exp(rxt)                   # 
#------------------------------------------#

Lam <- function(t1, x1, r1, k1, tau1){
    ( exp( t1*(r1*x1+k1) ) -1)/ ( exp( tau1*( r1*x1+k1 )  )-1 ) 
}

f <- function(t2, y2, u2, x2, r2, k2, tau2)
{ 
    Lam(t1=t2, x1=x2, r1=r2, k1=k2, tau1=tau2) - 
      u2 * Lam(t1=y2, x1=x2, r1=r2, k1=k2, tau1=tau2) 
}

# Lam(t1=0, x1=x, r1=r, k1=k, tau1=tau)-u* Lam(t1=y, x1=x, r1=r, k1=k, tau1=tau)

#f(t2=0, y2=y, u2 =u, x2=x, r2=r, k2=k, tau2=tau)

gendata<-function( N=200, a, b, r, k, tau=10 ){

    redata<- lwdata <-NULL
    for( i in 1:N){

        z <- rgamma(1, 4, 4) 
        x <- rbinom(1, 1, 0.5)
        #y <- runif(1, 0, 10)
        y <- min(rexp(1, 0.2), 10)
        m <- rpois(1, z* a * Lam(t1=y, x1=x, r1=r, k1=k, tau1=tau) *exp(b*x) ) 

        if (m>0){  	   	
            foo <- NULL
            for(mm in 1:m){ 
                u <- runif(1, 0, 1)
                foo <- c(foo, uniroot(f, lower=0, upper=y, y2=y, u2 =u, x2=x, r2=r, k2=k, tau2=tau)$root )
            }  	   	   
            foo<- sort(foo)

            foo1 <- c(0, foo)
            foo2 <- c(foo, y)
            foo3 <- c( rep(1, m), 0)
        } else { 
            foo <- foo1 <-foo3 <- 0 
            foo2 <- y }

        redata <- rbind( redata, cbind( id=i, t=foo, y=y, x=x, m=m) ) 
        lwdata <- rbind( lwdata, cbind( id=i, start=foo1, stop=foo2, event = foo3, x=x ) )
    } 

    list(re=redata, lw=lwdata) 

}



#========================================================#
#     calculate $\delta_{ijkl} * -(Xi-Xk)(tij-tkl)       # 
#                  assuming X(t) = X*t                   #       
#========================================================#   


R.cal <- function( datai, datak)
{
    yik <- min( datai$y[1], datak$y[1] )
    if ( any(datai$t <= yik) & any(datak$t <= yik) ){
        ti <- datai$t[ datai$t <= yik ]
        tk <- datak$t[ datak$t <= yik ]
        ( datai$x[1]- datak$x[1] ) * outer(ti, tk, "-") 
    }  else -99999 
} 


#========================================================# 
#                     calculate hik                      #
#                   assuming X(t) = X*t                  #        
#========================================================#         

hik <- function ( r,  R=R ) 
{
    h <- tapply( R$r, R$count,
                function(x, r0){
                    sum( x * exp( -r0*x )/(1 + exp(- r0*x) ) ) 
                }, r0=r)

    dh <- tapply( R$r, R$count,
                 function(x, r0){
                     sum( -x^2 * exp( -r0*x )/(1 + exp(-r0*x) )^2 ) 
                 }, r0=r) 
    cbind(  R[ cumsum( table(R$count) ), 2:3], 
          h= h,  dh=dh )
}



SStest <- function(id, t, y, x, weight="unit" )
{
    freq <- table(id)
    N <- length(freq)
    cfreq1 <- as.vector( cumsum(c(1,freq)) )
    cfreq2 <- as.vector( cumsum(freq) )
    Y.uniq<- tapply( y, id, mean)
    X.uniq<- tapply( x, id, mean)

    if (weight=="unit")  wt <- rep(1,N) else
        if (weight=="inverse mij") wt <- 1/table(id)
    #will not make a difference in deriving h if mi=0

    h.sh <- h.sz <- matrix(0, ncol=N, nrow=N)

    for (i in 2:N)
    {
        t1 <- t[cfreq1[i]:cfreq2[i]]
        for (j in 1:(i-1))
        {
            t2 <- t[cfreq1[j]:cfreq2[j]]
            y12 <- min( Y.uniq[i], Y.uniq[j])
            #assuming unit weight
            out <- sgnsum(t1=t1, t2=t2, y12=y12)
            h.sh[i,j] <- h.sh[j,i] <- wt[i]*wt[j]*sign(X.uniq[i]-X.uniq[j])*out$ss
            h.sz[i,j] <- h.sz[j,i] <- wt[i]*wt[j]*sign(X.uniq[i]-X.uniq[j])*out$diff
        }
    }
    ##Compute statistics
    w.sh <- sum(h.sh)/N/(N-1)
    w.sz <- sum(h.sz)/N/(N-1)

    ##Compute asymptotic variance of w.sh
    index <- 1:N
    Cov.sh <-Cov.sz <-0
    for( i in 1:N)
    {
        db.index <- expand.grid( index[-i],index[-i] )
        db.index <- db.index[ db.index[,1]!=db.index[,2],]
        Cov.sh <- Cov.sh +sum( h.sh[i,db.index[,1]] * h.sh[i,db.index[,2]])
        Cov.sz <- Cov.sz +sum( h.sz[i,db.index[,1]] * h.sz[i,db.index[,2]])
    }
    Cov.sh <- Cov.sh/N/(N-1)/(N-2)
    Cov.sz <- Cov.sz/N/(N-1)/(N-2)

    sigma.sh <-2* sqrt(Cov.sh/N)
    sigma.sz <-2* sqrt(Cov.sz/N)
    write.test(sh=w.sh, std.sh=sigma.sh, sz=w.sz, std.sz=sigma.sz)
    list(n=N,
         w.sh=w.sh, sigma.sh=sigma.sh,
         w.sz=w.sz, sigma.sz=sigma.sz,
         reject.sh= ( abs(w.sh)/sigma.sh>qnorm(0.975) ),
         reject.sz= ( abs(w.sz)/sigma.sz>qnorm(0.975) ) )
}



#   lw.est( x=lwdata$x; event=lwdata$event; start=lwdata$start; stop=lwdata$stop)

#  lw.est <- function(x, event, start, stop){
#  	
#  	   ut <- unique( stop[event==1] )      
#  	   
#  	   
#  	   
#  	   ss <- sum(x[event==1]) 
#  	   
#  	   U <- function(b, ut1=ut, x1=x; event1=event; start1=start; stop=stop; ss1=ss){
#  	   	
#  	   	       start
#  	   	}
#  	   
#  	    
#  	}

