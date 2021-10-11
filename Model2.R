################################################################################
# Model 2: Linear model to evaluate CpGs associated with change in PTSS in each cohort
################################################################################
library(lmerTest)
library(data.table)

# Scale pheno
pheno <- read.csv(" ")
pheno$PCL_SUM_01scaled <- (pheno$PCL_SUM-min(pheno$PCL_SUM))/(max(pheno$PCL_SUM)-min(pheno$PCL_SUM))

# Load beta matrix 
beta.norm<-fread(" ",data.table = F) # loading beta values, CpG X Sample
rownames(beta.norm)<-beta.norm$V1
beta.norm<-beta.norm[,-1]

# Define Variables
ageVar<-"AgeDiff" # Age variable = Age.post - Age.pre
cellTypes<-c("CD8T.Epic.diff","CD4T.Epic.diff","NK.Epic.diff","Bcell.Epic.diff","Mono.Epic.diff") # Cell type variables, difference in cell proportions between two time points
pcs<-c("Comp.2","Comp.3")  # Ancestry PC variables 
ptsdVar<-" " # PTSD variable = PTSS.post - PTSS.pre
studyID<-" " # E.g. "MRS","Prismo"
pheno$preID<-as.character(pheno$preID)
pheno$postID<-as.character(pheno$postID)

all(pheno$preID%in%colnames(beta.norm))  # Should be TRUE
all(pheno$postID%in%colnames(beta.norm)) # Should be TRUE

# Converting to M-values
range(beta.norm, na.rm=T)

# Changing beta values of 0 to 0.0001
if(min(beta.norm, na.rm=T)==0){
  beta.norm[which(beta.norm==0)]<-0.0001
}

# Changing beta values of 1 to 0.9999
if(max(beta.norm, na.rm=T)==1){
  beta.norm[which(beta.norm==1)]<-0.9999
}

range(beta.norm, na.rm=T)
sum(is.na(beta.norm))

# Convert to Mvalues using log2
beta.norm<-log2(beta.norm/(1-beta.norm)) # log transforming
range(beta.norm, na.rm=T)

formula<-as.formula(paste("outCpG~", 
                          paste(c("expCpG", ptsdVar, ageVar, cellTypes, pcs), collapse="+"), 
                          sep=""))
vars<-c("(Intercept)", "expCpG", ptsdVar, ageVar, cellTypes, pcs)
resultsBeta<-matrix(nrow=nrow(beta.norm), ncol=length(vars))
rownames(resultsBeta)<-rownames(beta.norm)
colnames(resultsBeta)<-vars
resultsSE<-resultsT<-resultsP<-resultsBeta
resultsDF<-matrix(nrow = nrow(beta.norm),ncol = 1)
rownames(resultsDF)<-rownames(beta.norm)
colnames(resultsDF)<-"df"

cpgs<-rownames(beta.norm)
errorProbes<-NULL

start<-proc.time()[3]
for(ii in 1:length(cpgs)){
  tempPheno<-pheno
  outCpG<-t(beta.norm[cpgs[ii], pheno$postID])
  expCpG<-t(beta.norm[cpgs[ii], pheno$preID])
  fit<-lm(formula, data=tempPheno)
  res<-coef(summary(fit))
  resultsBeta[ii,]<-res[, "Estimate"]
  resultsSE[ii,]<-res[, "Std. Error"]
  resultsT[ii,]<-res[, "t value"]
  resultsP[ii,]<-res[, "Pr(>|t|)"]
  resultsDF[ii,]<-fit$df.residual
  if(ii%%100==0){print(ii)}
}
end<-proc.time()[3]
end-start

beta<-as.data.frame(resultsBeta)
pval<-as.data.frame(resultsP)
se<-as.data.frame(resultsSE)
t<-as.data.frame(resultsT)

final<-as.data.frame(cbind(beta[ptsdVar],se[ptsdVar],t[ptsdVar],pval[ptsdVar],resultsDF))
rownames(final)<-rownames(beta)
colnames(final)<-c("BETA","SE","t","pval","df")

write.csv(final,file = paste0(studyID,"_Rutten_noob_allPreAsControls_PCL01scaled_covar_age_epicCellTypes_methPC.csv"),quote = F,row.names = F)




