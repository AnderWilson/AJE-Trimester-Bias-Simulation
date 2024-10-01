rm(list=ls())
gc()

library(ggplot2)
library(splines)
library(data.table)

# code with some helper functions to setup and analyze simulation
source("SimulationSetupCode.R")

# function to make plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



jackse <- function(x){
  n <- length(x)
  u <- rep(NA,n)
  for(i in 1:n){
    u[i] <- mean(x[-i])
  }
  return(sqrt(((n-1)/n)*sum((u-mean(u))^2)))
}


#--------------------------------------------------------------------------------------
# simulation

sderr <- 100

results <- data.table()
curves <- data.table()
for(simnum in 1:1000){
  message(simnum)
  set.seed(8657*simnum)
  res <- rnorm(n)
  
  for(fnum in 1:6){
    
    # simulate data
    xf <- linpred[[fnum]]
    y <- res*sderr  +xf
    
    
    # Fit joint TAE model and record results
    fit.together <- lm(y~xtri)
    fit.together.est <- c(fit.together$coef[-1],sum(fit.together$coef[-1])) 
    fit.together.se <- c(sqrt(diag(vcov(fit.together))[-1]),sqrt(sum(vcov(fit.together)[-1,-1]))) 
    results <- rbind(results,data.table(model="Joint TAE", trimester=trinames,est=fit.together.est, se=fit.together.se, df= fit.together$df.residual, simnum=simnum, f=fnum, dfspline=NA))
    curv <- c(rep(fit.together.est[1],12),rep(fit.together.est[2],14),rep(fit.together.est[3],11))
    curves <- rbind(curves,data.table(cbind(data.frame(model="Joint TAE",simnum=simnum, f=fnum, dfspline=NA),matrix(curv,1,37))))
    
    
    
    # fit separate TAE model and record results
    fit.separate.est <- c(lm(y~xtri[,1])$coef[-1],lm(y~xtri[,2])$coef[-1],lm(y~xtri[,3])$coef[-1],NA)
    fit.separate.se <- c(sqrt(diag(vcov(lm(y~xtri[,1])))[-1]),sqrt(diag(vcov(lm(y~xtri[,2])))[-1]),sqrt(diag(vcov(lm(y~xtri[,3])))[-1]),NA)
    results <- rbind(results,data.table(model="Separate TAE", trimester=trinames,est=fit.separate.est, se=fit.separate.se, df= n-2, simnum=simnum, f=fnum, dfspline=NA))
    curv <- c(rep(fit.separate.est[1],12),rep(fit.separate.est[2],14),rep(fit.separate.est[3],11))
    curves <- rbind(curves,data.table(cbind(data.frame(model="Separate TAE",simnum=simnum, f=fnum, dfspline=NA),matrix(curv,1,37))))
    
    # fine model fit statistics
    aic <- bic <- gcv <- NULL
    for(df in 3:10){
      mse <- mean((HI[[paste0("df",df)]]%*%y)^2)
      bic <- c(bic,n*log(mse)+log(n)*(1+df))
      aic <- c(aic,n*log(mse)+2*(1+df))
      gcv <- c(gcv,mse/tr[[paste0("df",df)]])
    }

    # choose best df based on one model fit statistic
    df <- c(3:20)[which.min(gcv)]
    
    # fit DLM and record results
    fit.ns <- lm(y~x%*%ns(1:37, df=df,intercept=TRUE))
    A <- lcomb%*%ns(1:37, df=df,intercept=TRUE)
    results <- rbind(results,data.table(model="DLM", trimester=trinames,est=drop(A%*%fit.ns$coef[-1]), se=sqrt(diag(A%*%vcov(fit.ns)[-1,-1]%*%t(A))), df= fit.ns$df.residual, simnum=simnum, f=fnum, dfspline=df))
    curv <- drop(ns(1:37, df=df,intercept=TRUE)%*%fit.ns$coef[-1])
    curves <- rbind(curves,data.table(cbind(data.frame(model="DLM",simnum=simnum, f=fnum, dfspline=NA),matrix(curv,1,37))))
    
  }
}

#save(results,curves,file="SimEstimates.Rdata")


#--------------------------------------------------------------------------------------
# summarize estimates curves
colnames(curves)[-c(1:4)] <- paste0("week",colnames(curves)[-c(1:4)])
curves <- data.table(curves)
setkeyv(curves,c("model","simnum","f"))
curves_fits <- melt(curves[,lapply(.SD, mean),by=c("model","f"),.SDcols=paste0("week",1:37) ], id.vars=c("model","f"))[,week:=as.numeric(substring(variable,5,6))]

# divide TAE estimates by number of weeks in trimester for making plots
# this divides the effect evently across weeks in the trimester.
curves_fits[model%in%c("Joint TAE","Separate TAE") & week%in%c(1:12), value:=value/12]
curves_fits[model%in%c("Joint TAE","Separate TAE") & week%in%c(13:26), value:=value/14]
curves_fits[model%in%c("Joint TAE","Separate TAE") & week%in%c(27:37), value:=value/11]

# format data generating mechanism data
fdat_dt <- data.table(fdat)[,list(f=GeneratingFunction, week=t, value=y)][,model:="True"][,variable:=paste0("week",week)][,f:=as.numeric(f)]

# format scenario name
curves_fits[,Scenario:=paste0("Scenario ",f)]
fdat_dt[,Scenario:=paste0("Scenario ",f)]


#--------------------------------------------------------------------------------------
# Final Publication Main Text Figure 2

# format scenario with letters
setkeyv(curves_fits, c("f","model","week"))
lets <- unique(curves_fits[,list(model,f)])
# strange letter arrangement below to change plot column order to match data anlaysis
lets[,panel:=LETTERS[1:nrow(lets)][rep(c(3,1,2),6)+(3*rep(0:5,each=3))]]
setkeyv(lets, c("f","model"))
curves_fits <- lets[curves_fits]

setkeyv(fdat_dt,  c("f","week"))
setkeyv(lets, c("f"))

fdat2 <- rbind(
  lets[model=="Joint TAE"][fdat_dt],
  lets[model=="Separate TAE"][fdat_dt],
  lets[model=="DLM"][fdat_dt])


fdat2[,Type:="Truth"]
curves_fits[,Type:="Estimated"]
plotdat <- rbind(fdat2[,list(model, f, panel, week, value, Type, Scenario)],
                 curves_fits[,list(model, f, panel, week, value, Type, Scenario)])


# make plots. 
# each panel is separate per journal request
plots <- list()

for(p in as.character(plotdat[,unique(panel)])){
fntsz <- 8
p_publish <- ggplot(data=plotdat[panel==p],aes(x=week, y=value, color=Type))
p_publish <- p_publish + ylab("Exposure effect") + xlab("Week")
p_publish <- p_publish + geom_line(size=1)
p_publish <- p_publish + theme(panel.background = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text  = element_text(size=fntsz, color="black"), 
                     axis.title = element_text(size=fntsz, color="black"),
                     legend.text = element_text(size=fntsz,color="black"),
                     legend.position="none", 
                     axis.line.x = element_line(color="black", size = 1/ggplot2:::.pt),
                     axis.line.y = element_line(color="black", size = 1/ggplot2:::.pt),
                     axis.ticks = element_line(colour = "grey20", size = 1/ggplot2:::.pt), 
                     axis.ticks.length = unit(fntsz/2, "pt"))
p_publish <- p_publish + ylim(-.32,1.22)  + scale_color_manual(values = c("black","darkgrey"))
p_publish <- p_publish + geom_hline(yintercept=0, lty=3, size = 1/ggplot2:::.pt)
p_publish <- p_publish + theme(plot.margin=unit(c(.1,0,0,.4), "cm"))
p_publish <- p_publish + theme(plot.title = element_text(hjust = -.38, size=fntsz,color="black")) + ggtitle(paste0(p,")"))

ggsave(filename=paste0("Figure2",p,".pdf"), plot=p_publish, path="FiguresAndTables",
       height=1.25, width=2)

ggsave(filename=paste0("Figure2",p,".eps"), plot=p_publish, path="FiguresAndTables",
       height=1.25, width=2, encoding="Greek")

plots[[p]] <- p_publish
}

# make combined version for convenience
pdf(file="FiguresAndTables/Figure2.pdf",height=1.25*6, width=6)
multiplot(plotlist = plots, layout = matrix(order(names(plots)),6,3,byrow=T))
dev.off()





#-------------------------------------------------------------------
# Bias in cumulative effect estimates for Web Table 1

# calculate true cumulative effect for each trimester
f <- as.numeric(substring(rownames(trisums),2,2))
trisums  <- data.table(trisums)
trisums[,f:=f]
tritrue <- melt(trisums, id="f", variable.name="trimester", value.name="truth")
setkeyv(tritrue,c("f","trimester"))

# summarize simulation results
results[,tot:=sum(est, na.rm=TRUE), by=c("model","f","simnum")]
results[is.na(est),est:=tot]
trismall <- results[,list(f,trimester,est,model)]
setkeyv(trismall,c("f","trimester","model"))

# merge simulation results and truth
cumulativetable <- trismall[tritrue]

# make bias variable
cumulativetable[, bias:=est-truth]

# results in long format
summarycum <- cumulativetable[,list(mean=mean(bias),bias=jackse(bias)),by=c("f","model","trimester")]

# make and save table in wide format
cumulative_table <- dcast(summarycum,f+model~trimester, value.var = c("mean","bias"))[, list(f,model,round(mean_Tri1,2),round(bias_Tri1,2),round(mean_Tri2,2),round(bias_Tri2,2),round(mean_Tri3,2),round(bias_Tri3,2),round(mean_All,2),round(bias_All,2))]
write.csv(cumulative_table,file="FiguresAndTables/EstimateCumulative.csv")
