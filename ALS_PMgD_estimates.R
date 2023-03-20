#### Calculate variant frequency in people affected by disease for ALS case studies

## Function for calculating PMgD from existing variant frequency estimates in familial and sporadic ALS and the 95% CI lower bound
PMgD_fromMetaanalysis<- function(PMgD_fam,
                                 PMgD_fam_95lower,
                                 PMgD_spor,
                                 PMgD_spor_95lower){
  
  #Determine PMgD as weighted sum reflecting ALS variant frequencies
  PMgD<- PMgD_fam*0.05+PMgD_spor*0.95
  
  #Determine 95% error and propagate the error summing in quadrature
  PMgD_fam_error <- PMgD_fam-PMgD_fam_95lower
  PMgD_spor_error <- PMgD_spor-PMgD_spor_95lower
  errBound <- sqrt((PMgD_fam_error*0.05)^2+(PMgD_spor_error*0.95)^2)
  
  #Write out the result with 95% CI
  writePretty<- formatC(
    c(PMgD,
      PMgD-errBound,
      PMgD+errBound),3,format="g")
  
  #Return nicely formatted
  return(paste0(writePretty[1], " (",writePretty[2],", ",writePretty[3],")"))
}


#SOD1 (all) estimate
PMgD_fromMetaanalysis(PMgD_fam=0.148,
                      PMgD_fam_95lower=0.115,
                      PMgD_spor=0.012,
                      PMgD_spor_95lower=0.007)
# "0.0188 (0.0138, 0.0238)"

#FUS (all) estimate
PMgD_fromMetaanalysis(PMgD_fam=0.028,
                      PMgD_fam_95lower=0.021,
                      PMgD_spor=0.003,
                      PMgD_spor_95lower=0.001)
#"0.00425 (0.00232, 0.00618)"

#C9orf72
PMgD_fromMetaanalysis(PMgD_fam=0.32,
                      PMgD_fam_95lower=0.28,
                      PMgD_spor=0.05,
                      PMgD_spor_95lower=0.04)
# "0.0635 (0.0538, 0.0732)"

















#SE of a proportion
se_p <- function(count,size){
  p <- count/size
  n <- size
  
  sqrt(p*(1-p)/n)
}



PMgD_fromCounts<- function(varFam,
                                 nFam,
                                 varSpor,
                                 nSpor){
  
  #Obtain Variant frequency estimates and the 95% error bounds in the familial and sporadic states
  PMgD_fam<- varFam/nFam
  PMgD_fam_error <- se_p(varFam,nFam)*1.96
  
  PMgD_spor <- varSpor/nSpor
  PMgD_spor_error <- se_p(varSpor,nSpor)*1.96
  #Obtain 95% error bounds
  
  #Estimate P(M|D) as weighted sum
  PMgD_total<- PMgD_fam*0.05+PMgD_spor*0.95
  
  #Propagate error
  errBound <- sqrt((PMgD_fam_error*0.05)^2+(PMgD_spor_error*0.95)^2)
  
  #Write out the result with 95% CI
  writePretty<- formatC(
    c(PMgD_total,
      PMgD_total-errBound,
      PMgD_total+errBound),3,format="g")
  paste0(writePretty[1], " (",writePretty[2],", ",writePretty[3],")")
}

#### FUS ClinVar
#Read in FUS per-variant counts
fusVar<- read.csv("~/OneDrive - King's College London/PhD/PhD project/Screening paper/git.local.genetScreening/FUS_clinvar_variants.csv")

#Run function
PMgD_fromCounts(varFam=sum(fusVar$Var_familial,na.rm = TRUE),
                nFam=median(fusVar$N_familial,na.rm = TRUE),
                varSpor=sum(fusVar$Var_sporadic,na.rm = TRUE),
                nSpor=median(fusVar$N_sporadic,na.rm = TRUE)
)
# "0.00251 (0.00125, 0.00377)"


#### SOD1 A5V

#Run function
PMgD_fromCounts(varFam=7,
                nFam=1125,
                varSpor=1,
                nSpor=4366
)
#"0.000529 (4.43e-05, 0.00101)"

