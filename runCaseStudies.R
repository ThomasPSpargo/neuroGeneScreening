######
### Apply genetic screening framework to estimate risk across case studies
### Written to the file output.csv; small additional calculations at the end of script
######

#Set wd to folder containing data and function
#setwd("path/to/repo")

#Run the function script
source("bayesScreening_func.R")

inputs <-read.csv("inputs.csv",header=T) #retrieve csv file
rownames(inputs) <- inputs[,1] #Assign rownames to dataset
inputs <- inputs[,-1] #rebuild object without column 1

outputs <- matrix(nrow=nrow(inputs),ncol=4) #Prepare object to store output
colnames(outputs) <- c("PDgT", "PDgT'","RR","PercentR") #Define column names
rownames(outputs) <- rownames(inputs) #Define row names by those of inputs

#Run the framework for all case study iterations
for(i in 1:nrow(inputs)){
  
  #Run the genetic screening
  runScreener<- bayesScreening(PD=inputs[i,1], PMgD=inputs[i,2],PDgM=inputs[i,3],TPR=inputs[i,4],TNR=inputs[i,5])
  
  #PDgT
  outputs[i,1]<- runScreener$PDgtestpos
  #PDgT'
  outputs[i,2]<- runScreener$PDgtestneg
  #RR
  outputs[i,3]<- outputs[i,1]/outputs[i,2]
  #% risk increase
  outputs[i,4]<- (outputs[i,3]-1)*100
}

outputs_sig<- signif(outputs, digits=3)
write.csv(outputs_sig,"outputs.csv") #Output the values of table 3 into a csv


######
### Run additional scenarios
######

#Relative disease risk for C9 secondary test relative to original negative
C9_second<- outputs['C9orf72 (RP valdate)' ,"PDgT"]/outputs['C9orf72',"PDgT'"]
#6.34828
PKU_second<- 0.889/1e-6
#889000

###


####Calculate the performance of tandem mass spec for screening PKU
PD = inputs["PAH (genetic)","PD"]
PTgD = 1
PNTgND =.9995
PT = (PTgD*PD)+((1-PNTgND)*(1-PD))
PDgT = PD*PTgD/PT
PDgT #Prob disease given pos test: 0.1666806

PDgNT = PD*(1-PTgD)/(1-PT) 
PDgNT #Prob disease given neg test: 0
####
