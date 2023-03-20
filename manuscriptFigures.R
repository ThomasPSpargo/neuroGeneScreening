######
### This script generates the figures for the manuscript
######


#Quick function to add panel letters to the top right of a panel
#https://seananderson.ca/2013/10/21/panel-letters/
add_label_legend <- function(pos = "topleft", label, ...) {
  legend(pos, label, bty = "n", ...)
}

#Set wd to allow running of the function
#setwd("path/to/repo")


#Run the function script
source("bayesScreening_func.R")
#Function  name: bayesScreening

###############################
#Name parameters used in modelled case studies

#Retrieve input data for case studies (represented in table 1 of ms)
inputs <-read.csv("inputs.csv",header=T) #retrieve csv file
rownames(inputs) <- inputs[,1] #Assign rownames to dataset
inputs <- inputs[,-1] #rebuild object without column 1

#P(D)
HD = inputs[grep("HTT",rownames(inputs))[1],"PD"]
ALS = inputs[grep("SOD1",rownames(inputs))[1],"PD"]
PKU = inputs[grep("PAH",rownames(inputs))[1],"PD"]

#Hypothetical P(Ds)
hyp_1e5 = as.numeric(format(1e-5,scientific=FALSE))
hyp_5e5 = as.numeric(format(5e-5,scientific=FALSE))

#P(M|D)
HTT = inputs["HTT (blind)","PMgD"]
SOD1 = inputs["SOD1 (all)","PMgD"]
A5V = inputs["SOD1 (A5V)","PMgD"]
FUSall = inputs["FUS (all)","PMgD"]
FUSpath = inputs["FUS (ClinVar)","PMgD"]
C9orf72 = inputs["C9orf72","PMgD"]
PAH = inputs["PAH (genetic)","PMgD"]


#Assign general parameters
#TEST Parameters
PointM = .9996
PointNM = .9995

IndelM=.9962
IndelNM=.9971

STRM=.99
STRNM =.90

CNVgridssM = .289
CNVgridssNM = .959

CNVdupM = .102
CNVdupNM = .9233

PerfectTest = 1


########################

#FIGURE 1
#Heatmap to display otherwise 3D surface
#X and Y axis are  TPR and TNR
#Z Axis is PDgtestpos
library(dplyr)
library(ggplot2)
library(gridExtra) #Used for arranging the multi-panel plot

#Give fixed values to the first 3 terms
PDgM=1
PMgD=1

#Define a sequence of values for TPR and TNR
y <- x <- seq(0.0,1,0.02)

####PANEL A
z <- matrix(nrow=length(x),ncol=length(y)) #Create empty Z matrix
rownames(z) <- colnames(z) <- x #define X and Y values associated with each value of Z

PD=.5 #Define PD in the first figure, 50%
#For loop to calculate matrix of z across values of x and y
for (i in 1:length(x)) {
  for(j in 1:length(y)) {
    z[i,j] <- bayesScreening(PD,PMgD,PDgM,x[i],y[j])$PDgtestpos
  }
}

zmelt <- reshape2::melt(z) #Reshape z for ggplot

#Create  plot from which a global legend can be extracted
fakep<- ggplot(zmelt,aes(x=Var1, y=Var2, fill=value))+
  geom_raster(interpolate = T)+
  labs(fill="P(D|T)")+
  scale_fill_distiller(type="seq", palette="PuOr", direction=-1,guide="colourbar", limits=c(0,1), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))

#Extract the heatmap colour-bar legend from the first plot (this legend is universal across plots)
g <- ggplotGrob(fakep)$grobs
colourbar_legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]

#Build real plot for panel A - save into object p
p<- ggplot(zmelt,aes(x=Var1, y=Var2, fill=value, z=value))+
  geom_raster(interpolate=TRUE) +
  labs(fill="P(D|T)", subtitle="(A)", caption="P(D) = 0.5, P(D|M) = 1, P(M|D) = 1")+
  theme(legend.position="none", panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        plot.subtitle=element_text(hjust=1), plot.caption=element_text(hjust=0.5,size=8))+
  scale_fill_distiller(type="seq", palette="PuOr", direction=-1,guide="colourbar")+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1),expand = F)+
  scale_x_continuous(breaks=seq(0,1,by =.2))+
  scale_y_continuous(breaks=seq(0,1,by =.2))+
  xlab("Sensitivity (P(T|M))")+
  ylab("Specificity (P(T'|M'))")

######PANEL B
z <- matrix(nrow=length(x),ncol=length(x)) #Create empty Z matrix
rownames(z) <- colnames(z) <- x #define X and Y values associated with each value of Z

PD=.01 #Define PD in the second figure, 1%
#For loop to calculate matrix of z across values of x and y
for (i in 1:length(x)) {
  for(j in 1:length(y)) {
    z[i,j] <- bayesScreening(PD,PMgD,PDgM,x[i],y[j])$PDgtestpos
  }
}
zmelt <- reshape2::melt(z) #Reshape z for ggplot

#Build plot for panel B - save into object q
q<- ggplot(zmelt,aes(x=Var1, y=Var2, fill=value, z=value))+
  geom_raster(interpolate=TRUE) +
  labs(fill="P(D|T)", subtitle="(B)", caption="P(D) = 0.01, P(D|M) = 1, P(M|D) = 1")+
  theme(legend.position="none", panel.border=element_rect(colour = "black", fill=NA, size=0.5),
        plot.subtitle=element_text(hjust=1), plot.caption=element_text(hjust=0.5,size=8))+
  scale_fill_distiller(type="seq", palette="PuOr", direction=-1,guide="colourbar")+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1),expand = F)+ #, clip="off")+ #clip off to allow annotation outside plot range
  scale_x_continuous(breaks=seq(0,1,by =.2))+
  scale_y_continuous(breaks=seq(0,1,by =.2))+
  xlab("Sensitivity (P(T|M))")+
  ylab("Specificity (P(T'|M'))")

#Write to file
svg(height=90/25.4, width=159.2/25.4, 
    pointsize=12, file="figure1heatmap.svg",family="arial")

#Arrange into 3 panel grid (panel 3 is global legend), allocating 4:4:1 space
gridExtra::grid.arrange(p,q,colourbar_legend, ncol=3,widths=c(4,4,1))

dev.off()


#Figure 2
#X=P(D), Y=P(D|T)
#Model across levels of analytical validity for the defined screening tools
#Could also make it P(M|T),
#png(height=120, width=100, pointsize=8, filename="figure2.png", res=300,units="mm")

svg(height=120/25.4, width=100/25.4, pointsize=10, file="figure2.svg",family="arial")

#Build a single plot panel:
opar<-par(mfrow=c(1,1))
xx=c(0,1) #Define for X scale
yy=c(0,1) #Define for Y scale
xlabel="Pre-test disease probability (P(D))"
ylabel1="P(D|T)" 
plot(xx,yy,ty="n",xlab=xlabel,ylab=" ",main=" "
     , sub="P(D|M) = 1, P(M|D) = 1", cex.sub=0.8) #Define plot, sub text for unvaried terms, no 'main' title

rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "floralwhite") #define plot colouring
mtext(ylabel1,side=2,line=2,col=1) #label Y axis

#Add h and v lines along plot
abline(h=.2,lty=1, lwd=1,col="mistyrose2")
abline(h=.4,lty=1, lwd=1,col="mistyrose2")
abline(h=.6,lty=1, lwd=1,col="mistyrose2")
abline(h=.8,lty=1, lwd=1,col="mistyrose2")

abline(v=.2,lty=1, lwd=1,col="mistyrose2")
abline(v=.4,lty=1, lwd=1,col="mistyrose2")
abline(v=.6,lty=1, lwd=1,col="mistyrose2")
abline(v=.8,lty=1, lwd=1,col="mistyrose2")

#Plot legend to acompany plotlines
legend(.55,.45,legend=c("SNV","Indel","STRE","CNV del.", "CNV dup."),
       col=c("forestgreen","red","navy","pink","purple"), lty=c(1,1,1,1,1),lwd=3, title="Variant type") 

#Construct necessary plot lines with curve function
#Specify fixed terms
PDgM=1 #Penetrance
PMgD=1 #Variant frequency

#Specify sensitivity and specificity for each tool and then plot curve model across P(D)
#Additional steps where plot points are added
TPR=IndelM
TNR=IndelNM
curve.model <- curve(bayesScreening(x,PMgD,PDgM,TPR,TNR)$PDgtestpos,from=0.000001,to=1,col="red",lty=1,lwd=2,add=T)

TPR=PointM
TNR=PointNM
curve.model <- curve(bayesScreening(x,PMgD,PDgM,TPR,TNR)$PDgtestpos,from=0.000001,to=1,col="forestgreen",lty=1,lwd=2,add=T)

TPR=STRM
TNR=STRNM
curve.model <- curve(bayesScreening(x,PMgD,PDgM,TPR,TNR)$PDgtestpos,from=0.000001,to=1,col="navy",lty=1,lwd=2,add=T)


TPR=CNVgridssM
TNR=CNVgridssNM
curve.model <- curve(bayesScreening(x,PMgD,PDgM,TPR,TNR)$PDgtestpos,from=0.000001,to=1,col="pink",lty=1,lwd=2,add=T)

TPR=CNVdupM
TNR=CNVdupNM
curve.model <- curve(bayesScreening(x,PMgD,PDgM,TPR,TNR)$PDgtestpos,from=0.000001,to=1,col="purple",lty=1,lwd=2,add=T)

dev.off()
#END run script for  FIGURE 2 











#Figure 3

svg(height=200/25.4, width=159.2/25.4, 
  pointsize=12, file="figure3.svg",family="arial")

#Build a plot panel:
opar<-par(mfrow=c(2,2)) #####TS  for 2x2 frame 
xx=c(0,1) #Adapt for X scale
yy=c(0,1) #Adapt for Y scale
xlabel="Penetrance (P(D|M))" #Model across penetrance
ylabel1="P(D|T)" #Y axis is P(D|T)

#Plot1 details - ALS SNVs across penetrance
plot(xx,yy,ty="n",xlab=xlabel,ylab=" ",main=" "
     , sub="P(D) = 0.0033, P(T|M) = 0.9996, P(T’|M’) = 0.9995", cex.sub=0.8)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "floralwhite") #grey94 #bisque
mtext(ylabel1,side=2,line=2,col=1)

#Lines for plot
abline(h=.2,lty=1, lwd=1,col="mistyrose2")
abline(h=.4,lty=1, lwd=1,col="mistyrose2")
abline(h=.6,lty=1, lwd=1,col="mistyrose2")
abline(h=.8,lty=1, lwd=1,col="mistyrose2")

abline(v=.2,lty=1, lwd=1,col="mistyrose2")
abline(v=.4,lty=1, lwd=1,col="mistyrose2")
abline(v=.6,lty=1, lwd=1,col="mistyrose2")
abline(v=.8,lty=1, lwd=1,col="mistyrose2")


#Legend with hypothetical variant
legend(0,1,legend=c("Hypothetical variant: 1","Hypothetical variant: 0.5","Hypothetical variant: 0.1", expression(paste(italic("SOD1"), ": 0.0188")), expression(paste(italic("SOD1"), " (A5V): 0.000529")), expression(paste(italic("FUS"), ": 0.0043")), expression(paste(italic("FUS"), " (ClinVar): 0.00251"))),
       col=c("black","pink","orange","forestgreen","green","mediumblue","lightblue"), lty=c(3,3,3,1,1,1,1),lwd=3, title="Case: P(M|D)", cex=.5)

add_label_legend("topright", "(A)")


#Define parameters for all examples:
PD=ALS
TPR=PointM
TNR=PointNM

#ALS-SOD1
PMgD=SOD1

#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="forestgreen",lty=1,lwd=2,add=T)

#Plot penetrance estimate on line
PDgM=inputs[grep("SOD1 \\(all)",rownames(inputs)),"PDgM"]

#Locating points to plot location of Penetrance for the gene simulation is based upon
# Approximation, with spline function interpolating curve for greater granularity
curve.interpol <- spline(curve.model, n=1000,  method = "natural")
values1 <- curve.interpol$x
values2 <- curve.interpol$y
xlocus <- which(abs(curve.interpol$x-PDgM)==min(abs(curve.interpol$x-PDgM))) #This is the closest to the other function
xid <- values1[xlocus] #Store in new xid object - later storage adds length on this object
yid <- values2[xlocus] #as above, for yid


#ALS-SOD1-A5V
PMgD=A5V

#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="green",lty=1,lwd=2,add=T)

#Plot penetrance estimate on line
PDgM=inputs[grep("SOD1 \\(A5V)",rownames(inputs)),"PDgM"]

#Locating points to plot location of Penetrance for the gene simulation is based upon
# Approximation, with spline function interpolating curve for greater granularity
curve.interpol <- spline(curve.model, n=1000,  method = "natural")
values1 <- curve.interpol$x
values2 <- curve.interpol$y
locus <- which(abs(curve.interpol$x-PDgM)==min(abs(curve.interpol$x-PDgM))) #This is the closest to the other function
xid[length(xid)+1] <- values1[locus]
yid[length(yid)+1] <- values2[locus]

#ALS-FUS-ALL
PMgD=FUSall

#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="mediumblue",lty=1,lwd=2,add=T)

PDgM=inputs[grep("FUS \\(all)",rownames(inputs)),"PDgM"]

#Locating points to plot location of Penetrance for the gene simulation is based upon
# Approximation, with spline function interpolating curve for greater granularity
curve.interpol <- spline(curve.model, n=1000,  method = "natural")
values1 <- curve.interpol$x
values2 <- curve.interpol$y
locus <- which(abs(curve.interpol$x-PDgM)==min(abs(curve.interpol$x-PDgM))) #This is the closest to the other function
xid[length(xid)+1] <- values1[locus]
yid[length(yid)+1] <- values2[locus]


#ALS-FUS-ClinVar
PMgD=FUSpath

#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="lightblue",lty=1,lwd=2,add=T)

PDgM=.591
PDgM=inputs[grep("FUS \\(ClinVar)",rownames(inputs)),"PDgM"]

#Locating points to plot location of Penetrance for the gene simulation is based upon
# Approximation, with spline function interpolating curve for greater granularity
curve.interpol <- spline(curve.model, n=1000,  method = "natural")
values1 <- curve.interpol$x
values2 <- curve.interpol$y
locus <- which(abs(curve.interpol$x-PDgM)==min(abs(curve.interpol$x-PDgM))) #This is the closest to the other function
xid[length(xid)+1] <- values1[locus]
yid[length(yid)+1] <- values2[locus]

#Hypothetical variants
#No point here because no estimate to be made
# # #Specify parameters:
PMgD = 1
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="black",lty=3,lwd=2,add=T)

PMgD = .5
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="pink",lty=3,lwd=2,add=T)

PMgD = .10
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="orange",lty=3,lwd=2,add=T)

points(xid, yid, pch =23, lwd = 2) #Plot point for penetrance

#Switch STR expansion performance
#Construct Plot 2:
plot(xx,yy,ty="n",xlab=xlabel,ylab=" ",main=" "
     , sub="P(D) = 0.0033, P(T|M) = 0.99, P(T’|M’) = 0.90", cex.sub=0.8)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "floralwhite") #grey94 #bisque
mtext(ylabel1,side=2,line=2,col=1)

#Lines for plot
abline(h=.2,lty=1, lwd=1,col="mistyrose2")
abline(h=.4,lty=1, lwd=1,col="mistyrose2")
abline(h=.6,lty=1, lwd=1,col="mistyrose2")
abline(h=.8,lty=1, lwd=1,col="mistyrose2")

abline(v=.2,lty=1, lwd=1,col="mistyrose2")
abline(v=.4,lty=1, lwd=1,col="mistyrose2")
abline(v=.6,lty=1, lwd=1,col="mistyrose2")
abline(v=.8,lty=1, lwd=1,col="mistyrose2")

#Legend with hypothetical variants
legend(0,1,legend=c("Hypothetical variant: 1","Hypothetical variant: 0.5","Hypothetical variant: 0.1",expression(paste(italic("C9orf72"),": 0.0635"))),
       col=c("black","pink","orange","red"), lty=c(3,3,3,1),lwd=3, title="Case: P(M|D)", cex=.5)

add_label_legend("topright", "(B)")


#Change analytic validity
TPR=STRM
TNR=STRNM

#ALS-C9orf72
PMgD=C9orf72

#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="red",lty=1,lwd=2,add=T)

PDgM=inputs[grep("C9orf72",rownames(inputs))[1],"PDgM"]

#Locating points to plot location of Penetrance for the gene simulation is based upon
# Approximation, with spline function interpolating curve for greater granularity
curve.interpol <- spline(curve.model, n=1000,  method = "natural")
values1 <- curve.interpol$x
values2 <- curve.interpol$y
xlocus <- which(abs(curve.interpol$x-PDgM)==min(abs(curve.interpol$x-PDgM))) #This is the closest to the other function
xid <- values1[xlocus]
yid <- values2[xlocus]


#Plot the points after all lines are plotted
#Hypothetical variants:
PMgD = 1
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="black",lty=3,lwd=2,add=T)

PMgD = .5
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="pink",lty=3,lwd=2,add=T)

PMgD = .10
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="orange",lty=3,lwd=2,add=T)


points(xid, yid, pch =23, lwd = 2) #Plot point for penetrance


#Construct Plot 3:

#Build a plot panel:
plot(xx,yy,ty="n",xlab=xlabel,ylab=" ",main=" "
     , sub="P(M|D) = 1, P(T|M) = 0.9996, P(T’|M’) = 0.9995", cex.sub=0.8)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "floralwhite") #grey94 #bisque
mtext(ylabel1,side=2,line=2,col=1)


#Lines for plot
abline(h=.2,lty=1, lwd=1,col="mistyrose2")
abline(h=.4,lty=1, lwd=1,col="mistyrose2")
abline(h=.6,lty=1, lwd=1,col="mistyrose2")
abline(h=.8,lty=1, lwd=1,col="mistyrose2")

abline(v=.2,lty=1, lwd=1,col="mistyrose2")
abline(v=.4,lty=1, lwd=1,col="mistyrose2")
abline(v=.6,lty=1, lwd=1,col="mistyrose2")
abline(v=.8,lty=1, lwd=1,col="mistyrose2")

#Legend
legend(0,1,legend=c("ALS: 0.0033", "HD: 0.00041", "PKU: 0.0001", "Hypothetical disease: 0.00001","Hypothetical disease: 0.00005"),
       col=c("black","forestgreen","pink","orange","mediumblue"), lty=c(1,1,1,3,3),lwd=3, title="Case: P(D)", cex=.5)

add_label_legend("topright", "(C)")

#Constant terms:
PMgD=1
TPR=PointM
TNR=PointNM

PD=ALS
#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="black",lty=1,lwd=2,add=T)

PD=HD
#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="forestgreen",lty=1,lwd=2,add=T)

PD=PKU
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="pink",lty=1,lwd=2,add=T)

PD=hyp_1e5
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="orange",lty=3,lwd=2,add=T)

PD=hyp_5e5
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="mediumblue",lty=3,lwd=2,add=T)




#Build the 4th plot panel:
plot(xx,yy,ty="n",xlab=xlabel,ylab=" ",main=" "
     , sub="P(M|D) = 1, P(T|M) = 0.99, P(T’|M’) = 0.90", cex.sub=0.8)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "floralwhite") #grey94 #bisque
mtext(ylabel1,side=2,line=2,col=1)

#Lines for plot
abline(h=.2,lty=1, lwd=1,col="mistyrose2")
abline(h=.4,lty=1, lwd=1,col="mistyrose2")
abline(h=.6,lty=1, lwd=1,col="mistyrose2")
abline(h=.8,lty=1, lwd=1,col="mistyrose2")

abline(v=.2,lty=1, lwd=1,col="mistyrose2")
abline(v=.4,lty=1, lwd=1,col="mistyrose2")
abline(v=.6,lty=1, lwd=1,col="mistyrose2")
abline(v=.8,lty=1, lwd=1,col="mistyrose2")

#Legend
legend(0,1,legend=c("ALS: 0.0033", "HD: 0.00041", "PKU: 0.0001", "Hypothetical disease: 0.00001","Hypothetical disease: 0.00005"),
       col=c("black","forestgreen","pink","orange","mediumblue"), lty=c(1,1,1,3,3),lwd=3, title="Case: P(D)", cex=.5)


add_label_legend("topright", "(D)")


#Update to STR analytical validity
TPR=STRM
TNR=STRNM


#Lines
PD=ALS
#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="black",lty=1,lwd=2,add=T)

PD=HD
#Plot estimates at specified mutation frequency
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="forestgreen",lty=1,lwd=2,add=T)

PD=PKU
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="pink",lty=1,lwd=2,add=T)

PD=hyp_1e5
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="orange",lty=3,lwd=2,add=T)

PD=hyp_5e5
#Add plot line for P(D|T) across penetrance:
curve.model <- curve(bayesScreening(PD,PMgD,x,TPR,TNR)$PDgtestpos,from=0.01,to=1,col="mediumblue",lty=3,lwd=2,add=T)

#End plot 4 and the final panel

dev.off()

