#ex 1 from vignette
set.seed( 123 )
pData <- data.frame(
  id = rep( paste( "F", 1:15, sep = "_" ), each = 4 ),
  time = rep( 1981:1984, 15 ) )
pData$mu <- rep( rnorm( 15 ), each = 4 )
pData$x1 <- rnorm( 60 )
pData$x2 <- runif( 60 )
pData$ys <- -1 + pData$mu + 2 * pData$x1 + 3 * pData$x2 + rnorm( 60 )
pData$y <- ifelse( pData$ys > 0, pData$ys, 0 )
library( plm )
pData <- pdata.frame( pData, c( "id", "time" ) )

panelResult <- censReg( y ~ x1 + x2, data = pData, method = "BHHH",
                        left=0, right=Inf)
summary( panelResult )




## ex 2 
# https://stats.stackexchange.com/questions/544652/how-to-write-a-random-effect-in-tobit-models-in-r-using-the-censreg-package/582354#582354

x=rep(c(0,2,4,6,8),20)
beta0=rep(rnorm(20,0,1),each=5)
ind=rep(1:20,each=5)
beta=2
y=x*beta + beta0 +rnorm(100,0,0.5)
data=data.frame(x=x,y=y,ind=ind) %>% mutate(y_cens=ifelse(y>12,12,y))


pData <- pdata.frame(data, c( "ind", "x" ) )
pData$x1 <- as.numeric(pData$x)
fit <- censReg(y_cens ~ x1, data = pData, method = "BHHH", left=-Inf, right=12)
summary(fit)


