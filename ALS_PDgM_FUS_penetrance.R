######
# This script calculates penetrance estimates for the FUS gene screening scenarios
# adpenetrance is documented at: https://github.com/thomaspspargo/adpenetrance
######

######
### Load functions
######
# Run the adpenetrance function script from GitHub
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/adpenetrance_function.R")
source("https://raw.githubusercontent.com/ThomasPSpargo/adpenetrance/master/getResidualRisk.R")

# SE of a proportion
se_p <- function(count,size){
  p <- count/size
  n <- size
  
  sqrt(p*(1-p)/n)
}

#---------------------------------------------------------------#
#       Calculate penetrance For FUS (all) and FUS (ClinVar)    #
#---------------------------------------------------------------#

#Universal parameters
N=1.543
PF=0.05
PA=1/300

#FUS (ALL) -
MF=0.028
MS=0.003
adpenetrance(N=N,
             MF=MF, MF_SE=(MF-0.021)/1.96,
             MS=MS, MS_SE=(MS-0.001)/1.96,
             PF=PF,
             useG= getResidualRisk(PA=PA, MF = MF, MS = MS,MU=0, PF = PF)
)

#P(fam) = 0.3294 (0.1721, 0.4867)
#f      = 0.579 (0.291, 0.884)

#FUS (ClinVar) -
#SE derived per FUS_pathogenic_error_estimation file in parameters folder.
MF=17/1108
MS=8/4366
adpenetrance(N=N,
             MF=MF, MF_SE = se_p(17,1108),
             MS=MS, MS_SE = se_p(8, 4366),
             PF=PF,
             useG= getResidualRisk(PA=PA, MF = MF, MS = MS,MU=0, PF = PF)
)

#P(fam) = 0.3059 (0.1280, 0.4838)
#f      = 0.536 (0.211, 0.877)
