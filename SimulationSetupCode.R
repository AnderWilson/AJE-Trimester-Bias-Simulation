

#--------------------------------------------------------------------------------------
#create functions to generate data from
f1 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(0,10),rep(1,4)))
f2 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(0,6),rep(1,2),rep(0,6)))
f3 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(0,7),rep(1,5),rep(0,2)))
f4 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(0,7),rep(1,7)))
f5 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(.5,14)))
f6 <- drop(bs(1:37,knots=seq(1,37,length=12)[-c(1,12)],intercept=TRUE)%*%c(rep(0,14)))
fdat <- rbind(
  data.frame(GeneratingFunction=1,t=1:37,y=f1),
  data.frame(GeneratingFunction=2,t=1:37,y=f2),
  data.frame(GeneratingFunction=3,t=1:37,y=f3),
  data.frame(GeneratingFunction=4,t=1:37,y=f4),
  data.frame(GeneratingFunction=5,t=1:37,y=f5),
  data.frame(GeneratingFunction=6,t=1:37,y=f6))
fdat$GeneratingFunction <- as.factor(fdat$GeneratingFunction)
fdat$Scenario <- paste("Scenario",fdat$GeneratingFunction)


#--------------------------------------------------------------------------------------
#find trimester averages
lcomb <- matrix(0,4,37)
lcomb[1,1:12] <- 1
lcomb[2,13:26] <- 1
lcomb[3,27:37] <- 1
lcomb[4,] <- 1

trisums <- t(lcomb%*%cbind(f1,f2,f3,f4,f5,f6))
colnames(trisums) <- c("Tri1","Tri2","Tri3","All")


#--------------------------------------------------------------------------------------
# some exposure data

x <- read.csv("exposure_data.csv")
x <- as.matrix(x)
xtri <- cbind(rowMeans(x[,1:12]),
              rowMeans(x[,13:26]),
              rowMeans(x[,27:37]))


#--------------------------------------------------------------------------------------
#things to save computation in the loop
n <- nrow(x)
linpred <- list(x%*%f1,x%*%f2,x%*%f3,x%*%f4,x%*%f5,x%*%f6)

trinames <- c("Tri1","Tri2","Tri3","All")

#this is the H matrix for each DF to be tested
H <- HI <- tr <- list()
for(df in 2:10){
  D <- cbind(1,x%*%ns(1:37, df=df,intercept=TRUE))  
  H[[paste0("df",df)]] <- D%*%solve(t(D)%*%D)%*%t(D)
  HI[[paste0("df",df)]] <- H[[paste0("df",df)]]-diag(n)
  tr[[paste0("df",df)]] <- (1-sum(diag(H[[paste0("df",df)]]))/n)^2 
}

