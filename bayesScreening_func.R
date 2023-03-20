######
### Function to calculate disease risk following positive or negative genetic screening result
######

bayesScreening<-function(PD,PMgD,PDgM,TPR,TNR){
  ########### #Define FNR and FPR for simplicity
  FNR= 1-TPR        # P(T-|M+)
  FPR= 1-TNR        # P(T+|M-)
  Rel= (TPR+TNR)/(TPR+TNR+FPR+FNR) #Reliability of test = accuracy =(TPR+TNR)/2
  ##########
  
  #Draft Eqn 1
  PM=PD*PMgD/PDgM #Total Probability of variant
  
  #Draft Eqn 2
  PDgNM = PD*(1-PMgD)/(1-PM)  # Probability of disease given no variant
  
  #Draft Eqn 3
  Ptestpos=TPR*PM+FPR*(1-PM)  # Probability of positive test
  
  ## Calculating PD|T
  #Draft Eqn 4
  PMgtestpos=PM*TPR/Ptestpos  # Probability of variant given positive test
  
  #Draft Eqn 6
  PDgtestpos= PDgM*PMgtestpos + PDgNM*(1-PMgtestpos) # Probability of disease given positive test
  
  
  ## Calculating PD|T-
  #Draft Eqn 5
  PMgtestneg=PM*FNR/(1-Ptestpos)  # Probability of variant given negative test
  
  #Draft Eqn 7
  PDgtestneg=PDgM*PMgtestneg + PDgNM*(1-PMgtestneg) # Probability of disease given negative test
  
  return(list(Rel=Rel, Ptestpos=Ptestpos,PMgtestpos=PMgtestpos,PDgtestpos=PDgtestpos,
              PDgtestneg=PDgtestneg,PMgtestneg=PMgtestneg, PD=PD))
}
