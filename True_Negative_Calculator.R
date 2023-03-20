#### Calculate sensitivity for the NGS genotyping cools


#Calculates True/False Positive and True/False negative rates based on Precison, recall and F1 scores (N=2)
posnegcalc <- function(prec, recall){
  
  #FP=B/A-B

  TP = recall
  FP=(TP/prec)-TP #Alternative formulation: FP = (TP-prec*TP)/prec
  FN = 1 -TP
  TN = 1- FP
  Mtot = TP + FN
  NMtot = FP + TN 
  
  #Old formulation: function(N, prec, recall, F1)
  #TP = recall                          #By definition
  #FP = (TP-prec*TP)/prec               #PPV formula
  #FN = (2*TP - (F1*2*TP + F1*FP))/F1   #F1 score formula
  #TN = N - (TP + FP + FN)              

  confusion_matrix<- c(TP, FP, FN, TN, Mtot, NMtot)
  confusion_matrix<- matrix(confusion_matrix, nrow = 2, ncol =3)
  colnames(confusion_matrix) <- c("TestPos", "TestNeg", "ConditionTotal")
  rownames(confusion_matrix) <- c("ConditionPos", "ConditionNeg")

  confusion_matrix
}

#----------------#
#   Get params   #
#----------------#

#Dragen Pipeline (SNVs)
posnegcalc(prec=0.9995, recall=0.9996)

#Dragen Pipeline (Indels)
posnegcalc(prec=.9971, recall=.9962)

#ExpansionHunter
posnegcalc(prec=0.91, recall=.99)

#GRIDSS
posnegcalc(prec=.876, recall=.289)

#Wham
posnegcalc(prec=.571, recall=.102)
