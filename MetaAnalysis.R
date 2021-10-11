########################################################################################
# PGC Meta-Analysis: CpG Sites in All Studies
########################################################################################

########################################################################################
# Step 1: Load and subset data
########################################################################################
# Step 1A: Load summary statistics from each cohort
starrs.results <- read.csv(" ", row.names = 1)
MRS.results <- read.csv(" ", row.names = 1)
PRISMO.results <- read.csv(" ",row.names = 1)

# Step 1B: Determine the number of probes available in all studies
starrs.sites <- rownames(starrs.results) 
MRS.sites <- rownames(MRS.results) 
PRISMO.sites <- rownames(PRISMO.results) 

sites<-append(starrs.sites, MRS.sites)
sites <- append(sites,PRISMO.sites)
sites<-unique(sites)
length(sites) 

all <- intersect(starrs.sites, MRS.sites)
all <- intersect(all,PRISMO.sites)

sum(is.na(match(all, rownames(starrs.results)))) # these should all be 0 
sum(is.na(match(all, rownames(MRS.results))))
sum(is.na(match(all, rownames(PRISMO.results))))
length(all)  

rm(starrs.sites,MRS.sites,PRISMO.sites)

########################################################################################
# Step 2: Data Prep
########################################################################################
# Step 2A: Subset all sites
studies <- c("starrs","MRS","PRISMO")

for(ii in 1:length(studies)){
  results<-get(paste(studies[ii], ".results", sep=""))
  results<-results[all, ]
  assign(paste(studies[ii], ".results", sep=""), results)
  rm(results)
}

table(rownames(MRS.results)==rownames(starrs.results))
table(rownames(MRS.results)==rownames(PRISMO.results))

# Step 2B: Calculate one-sided p-values for each CpG site's t-statistic
starrs.oneSided <- pt(starrs.results$t, df = starrs.results$df)
MRS.oneSided <- pt(MRS.results$t, df = MRS.results$df)
PRISMO.oneSided <- pt(PRISMO.results$t, df = PRISMO.results$df)

########################################################################################
# Step 3: Meta-Analysis of Sites in All Studies
########################################################################################

# beta coefficients
betas<-matrix(nrow=nrow(MRS.results), ncol=length(studies))
colnames(betas)<-studies
rownames(betas)<-rownames(MRS.results)

for(ii in 1:ncol(betas)){
  temp<-get(paste(studies[ii], ".results", sep=""))
  betas[, studies[ii]]<-temp[, "BETA"] # beta coefficient
  rm(temp)
}

# t-statistics from 1 sided p-values
tstats<-matrix(nrow=nrow(MRS.results), ncol=length(studies))
colnames(tstats)<-studies
rownames(tstats)<-rownames(MRS.results)

for(ii in 1:ncol(tstats)){
  temp<-get(paste(studies[ii], ".oneSided", sep=""))
  tstats[, studies[ii]]<-qnorm(as.numeric(temp)) # t-statistic from 1-sided p-value
  rm(temp)
}

# SEs
ses<-betas/tstats
bweights<-1/(ses^2)

# Calculate inverse-variance beta coefficient and variance
betaC<-apply(betas*bweights, 1, sum)/apply(bweights, 1, sum)
variance<-1/apply(bweights, 1, sum)
rm(betas, ses, bweights)

# Weights (df+covariates+1)
weights<-matrix(nrow=nrow(MRS.results), ncol=length(studies))
colnames(weights)<-studies
rownames(weights)<-rownames(MRS.results)

weights[,"starrs"] <- starrs.results$df + 1 + 8
weights[,"MRS"] <- MRS.results$df + 1 + 8
weights[,"PRISMO"] <- PRISMO.results$df + 1 + 8

for(ii in 1:nrow(weights)){
  weights[ii, ]<-sqrt(weights[ii, ]/sum(weights[ii, ]))
}

# Calculated Z-score and p-value
w.tstats<-weights*tstats
z.combined<-apply(w.tstats, 1, sum)
p<-2*(1-pnorm(abs(z.combined)))

all(names(p)==names(z.combined))
all(names(p)==names(betaC))
all(names(p)==names(variance))
res<-as.data.frame(cbind(z.combined, p, betaC, variance))

write.csv(res, file = " ",row.names = F)


