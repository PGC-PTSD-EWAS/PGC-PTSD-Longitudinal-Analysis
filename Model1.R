################################################################################
# Model 1: RE model to evaluate CpGs associated with PTSS in each cohort
################################################################################

library(lme4)
library(lmerTest)
library(data.table)

# Load methylation data
beta.norm<-fread(" ",data.table=F) # beta matrix. CpG x Participants
rownames(beta.norm)<-beta.norm$V1
beta.norm<-beta.norm[,-1]

# Load phenotype file
pheno<-read.csv(" ",stringsAsFactors=FALSE,header=TRUE, row.names = 1)

# Define Variables
ptsdVariable<-" " # 0 -1 scaled
covariates<-c("age","CD8T.Epic","CD4T.Epic","NK.Epic","Bcell.Epic","Mono.Epic","Comp.2","Comp.3") # E.g., c("Age", "Gender", "CD8T", ....)
idVariable<-"studyid" # Participant ID for the Random Effects. E.g., "IID"
studyID<-" " # E.g. "MRS","Prismo"

beta.norm<-beta.norm[,names(beta.norm)%in%row.names(pheno)]
pheno<-pheno[row.names(pheno)%in%names(beta.norm),]
beta.norm<-beta.norm[,order(names(beta.norm))]
pheno<-pheno[order(row.names(pheno)),]
table(colnames(beta.norm)==rownames(pheno))
range(pheno[,ptsdVariable])

# Remove NAs
naindex<-pheno[,c(ptsdVariable,covariates,idVariable)]
naindex<-complete.cases(naindex) 
beta.norm<-beta.norm[,naindex]
pheno<-pheno[naindex,]
table(rownames(pheno)==colnames(beta.norm))

# Converting to M-values
range(beta.norm, na.rm=T)

# Changing beta values of 0 to 0.0001
if(min(beta.norm, na.rm=T)==0){
  beta.norm[which(beta.norm==0)]<-0.0001
}

# Changing beta values of 1 to 0.9999
if(max(beta.norm, na.rm=T)==1.000000e+00){
  beta.norm[which(beta.norm==1.000000e+00)]<-0.9999
}

range(beta.norm, na.rm=T)
sum(is.na(beta.norm))

# Convert to Mvalues using log2
beta.norm<-log2(beta.norm/(1-beta.norm)) # log transforming
sum(beta.norm=="-Inf", na.rm=T) # Should be 0
sum(is.na(beta.norm)) # should be same as the number of missing in the beta.norm matrix
range(beta.norm, na.rm=T)

# Configuring phenotypes
pheno$meth<-NA
pheno[,idVariable]<-factor(pheno[,idVariable])

# Results objects
resultsBeta<-matrix(nrow=nrow(beta.norm), ncol=2+length(covariates))
rownames(resultsBeta)<-rownames(beta.norm)
colnames(resultsBeta)<-c("(Intercept)", ptsdVariable, covariates)
resultsSE<-resultsT<-resultsP<-resultsDF<-resultsBeta
errorProbes<-NULL
warningProbes<-NULL

assign("last.warning", NULL, envir = baseenv())

formula<-as.formula(paste("meth~",ptsdVariable, "+", 
                          paste(covariates, collapse="+"), 
                          "+(1|", idVariable, ")", sep=""))
# Running analyses
start<-proc.time()[3]
for(ii in 1:nrow(beta.norm)){
  pheno$meth<-t(beta.norm[ii, rownames(pheno)])
  fit<-try(lmer(formula, data=pheno), silent=F)
  if(class(fit)!="try-error" & is.null(warnings())){
    res<-coef(summary(fit))
    resultsBeta[ii,]<-res[, "Estimate"]
    resultsSE[ii,]<-res[, "Std. Error"]
    resultsT[ii,]<-res[, "t value"]
    resultsP[ii,]<-res[, "Pr(>|t|)"]
    resultsDF[ii,]<-res[, "df"]
  }
  
  if(class(fit)=="try-error"){
    errorProbes<-append(errorProbes, rownames(beta.norm)[ii])
  }
  
  if(!is.null(warnings())){
    warningProbes<-append(warningProbes, rownames(beta.norm)[ii])
    assign("last.warning", NULL, envir = baseenv())
  }
  
  if(ii%%10==0){print(ii)}
  pheno$meth<-NA
}
end<-proc.time()[3]
end-start

rm(res, fit)

# For main Pheno
final<-as.data.frame(cbind(resultsBeta[,ptsdVariable],resultsSE[,ptsdVariable],resultsT[,ptsdVariable],resultsP[,ptsdVariable],resultsDF[,ptsdVariable]))
rownames(final)<-rownames(beta.norm)
colnames(final)<-c("BETA","SE","t","pval","df")

write.csv(final,file = paste0(studyID,"_RE_PCL01scaled_covar_age_epicCellTypes_methPC.csv"),quote = F,row.names = F)


