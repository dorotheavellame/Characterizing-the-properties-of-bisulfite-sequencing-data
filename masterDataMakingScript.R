#######################################################
### Scripts used to compile data in RRBS data paper ###
### Dorothea Seiler Vellame ###########################
### ds420@exeter.ac.uk ################################
#######################################################

### Data needed #######################################

# rTg4510MatrixUnfiltered.Rdata
# J20MatrixUnfiltered.Rdata
# Normalised_Data_Sesame.rdat
# rTg4510_phenotype_RRBS.csv
# J20_phenotype_RRBS.csv
# HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.1.csv

### Data made #########################################

# RRBSmergedThatCorrelate.Rdata
# mappingFull.Rdata
# arrayRRBSComparisonData.Rdata
# RDThresholdOutput.Rdata
# longDat.Rdata
# benchmarkingDataAll.Rdata


## merge RRBS data to make one big data set
## read in and correlate the seperate data 
load("/mnt/data1/Thea/IsabelSamples/data/rTg4510MatrixUnfiltered.Rdata")
rtg = RRBSmatrix
load("/mnt/data1/Thea/IsabelSamples/data/J20MatrixUnfiltered.Rdata")
j20 = RRBSmatrix

rm(RRBSmatrix)

## get overlapping CpGs (merge by col 1)
all = merge(rtg, j20, by = "chr_start")

## reorder all columns
allOrdered = all[, c(1, 126:188, 2:63, 189:251, 64:125)]

## save all as it doesn't have correlation issues and rerun
save(allOrdered, file = "/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")


load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
load("/mnt/data1/Dementia_Mouse/Array/Normalised_Data_Sesame.rdat")

rtgPheno = read.csv("/mnt/data1/Thea/IsabelSamples/data/rTg4510_phenotype_RRBS.csv")
j20Pheno = read.csv("/mnt/data1/Thea/IsabelSamples/data/J20_phenotype_RRBS.csv")

## rename data for ease
RRBS = allOrdered
phenoRRBS = rbind.data.frame(j20Pheno[,colnames(j20Pheno) %in% colnames(rtgPheno)], rtgPheno[,colnames(rtgPheno) %in% colnames(j20Pheno)])
array = Normalised_Sesame_Betas
phenoArray = QCmetrics

rm(allOrdered, j20Pheno, rtgPheno, Normalised_Sesame_Betas, QCmetrics)

## keep only the sites in both the array and RRBS
mappings = read.csv("/mnt/data1/Dementia_Mouse/MammalianArrayNormalizationTools/mappabilityData/HorvathMammalChip_SpeciesGenomeCoordsUnique_v1.1.csv")
mapping = mappings[,c("probeID", "MusMusculus")]

## change format of MusMusculus column to match RRBS data
x = strsplit(as.character(mapping[,2]), split = ":")

y = c()
for (i in 1:length(x)){
  if (is.na(x[[i]][1])){
    y = c(y, NA)
  }else{y = c(y, paste("chr", x[[i]][1], "_", x[[i]][2], sep = ""))
  }
}
mapping[,3] = y

rm(x, y)

## remove all those from mapping that don't have an RRBS value
mapping = mapping[!is.na(mapping[,3]),]


## create boolean columns in mapping for those CPGs in array, RRBS and both
mapping[,4] = FALSE
mapping[,5] = FALSE
mapping[which(mapping$probeID %in% rownames(array)), 4] = TRUE
mapping[which(mapping$V3 %in% RRBS$chr_start), 5] = TRUE
mapping[,6] = mapping[,4] & mapping[,5]

mapping = mapping[mapping[,6], c(1,3)]

colnames(mapping) = c("arrayID", "RRBSID")

## remove non unique CpGs
# remove non unique CpGs
if(sum(duplicated(mapping$RRBSID)) > 0){
  mapping = mapping[-which(duplicated(mapping$RRBSID)),]
}
save(mapping, file = "/mnt/data1/Thea/Simulations/data/ValidationData/mappingFull.Rdata")

## match array and RRBS order to mapping
arrayindex = match(mapping$arrayID, rownames(array))
arrayindex = arrayindex[!is.na(match(mapping$arrayID, rownames(array)))]
arraysub = array[arrayindex,]
RRBSindex = match(mapping$RRBSID, RRBS$chr_start)
RRBSindex = RRBSindex[!is.na(match(mapping$RRBSID, RRBS$chr_start))]
RRBSsub = RRBS[RRBSindex,]

## check that orders match mapping
sum(mapping$arrayID == rownames(arraysub)) == 3352
sum(mapping$RRBSID == RRBSsub$chr_start) == 3352


## check which samples overlap, subset and merge phenotype files
## 80 samples
phenoArraysub = phenoArray[1:80,]
x = strsplit(phenoArray$ExternalSampleID, split = "_")
y = c()
for (i in 1:length(x)){
  y = c(y, x[[i]][2])
}

colnames(arraysub) = y

phenoArraysub  = cbind.data.frame(phenoArraysub, y[1:80])
rm(x, y)

pheno = merge(phenoArraysub, phenoRRBS, by.x = "y[1:80]", by.y = "Sample_ID")

## remove samples that don't have corresponding RRBS data
arraysub = arraysub[,which(colnames(arraysub) %in% pheno$y)]


## split RRBS data and remove samples not in array
meth = cov = matrix(nrow = nrow(RRBSsub), ncol = 0)
name = as.character(pheno$y)

for (i in 1:length(name)){
  ## take name[i] column name and put first in meth
  meth = cbind(meth, RRBSsub[,grepl(name[i], colnames(RRBSsub))][,1])
  ## put second in cov
  cov = cbind(cov, RRBSsub[,grepl(name[i], colnames(RRBSsub))][,2])
}
RRBSmeth = apply(meth, 2, as.numeric)/100
RRBScov = apply(cov, 2, as.numeric)

RRBSsub = cbind.data.frame(RRBSsub$chr_start, RRBSmeth, RRBScov)
colnames(RRBSsub) = c("chr_start", paste(name, "_m"), paste(name, "_cov"))
colnames(RRBSmeth) = colnames(RRBScov) = name
rownames(RRBScov) = rownames(RRBSmeth) = RRBSsub[,1]

## save data
save(RRBScov, RRBSmeth, arraysub, pheno, RRBSsub, file = "/mnt/data1/Thea/Simulations/data/ValidationData/arrayRRBSComparisonData.Rdata")



load("/mnt/data1/Thea/Simulations/data/ValidationData/arrayRRBSComparisonData.Rdata")
covFilter = 25 # arbitrary choice to show effect of some amount of filtering

## melt the data and form one long data frame
arraylong = arraysub[,1]
RRBSmethlong = RRBSmeth[,1]
Coverage = RRBScov[,1]

for (i in 2:ncol(arraysub)){
  arraylong = c(arraylong, arraysub[,i])
  RRBSmethlong = c(RRBSmethlong, RRBSmeth[,i])
  Coverage = c(Coverage, RRBScov[,i])
}


longDat = data.frame(as.numeric(arraylong), as.numeric(RRBSmethlong), as.numeric(Coverage))
colnames(longDat) = c("arraylong", "RRBSmethlong", "Coverage")

## add column for samplename and genotype from pheno and model
longDat$Sample = rep(colnames(arraysub), each = nrow(arraysub))
longDat$Genotype = rep(pheno$Genotype.x, each = nrow(arraysub))
longDat$Model = rep(pheno$AD_model, each = nrow(arraysub))
rm(RRBScov, RRBSmeth, arraysub, RRBSsub)

## remove all rows with NA in RRBS data
longDat = longDat[!is.na(longDat$RRBSmethlong),]

## remove CpGs with coverage above 125 as too large is probably technical noise
longDat = longDat[which(longDat$Coverage < 125),]

## remove those with cov < covFilter
longDatFilt = longDat[which(longDat$Coverage > covFilter),]

save(longDat, longDatFilt, file = "/mnt/data1/Thea/Simulations/data/ValidationData/longDat.Rdata")



## make data for time vs accuracy of tool analysis
load("/mnt/data1/Thea/Simulations/data/ValidationData/RRBSmergedThatCorrelate.Rdata")
RRBSMatrix = allOrdered
covCols = 127:251
source("/mnt/data1/Thea/Simulations/RScripts/POWEREDBiSeqFunction.R")

nSamples = (ncol(RRBSMatrix)-1)/2

nSamplesPerGroup1 = floor(nSamples/2)
nSamplesPerGroup2 = ceiling(nSamples/2)
pheno = as.factor(c(rep("case", each = nSamplesPerGroup1),rep("control", each = nSamplesPerGroup2)))
rd = 20
nSampleNeeded = 2

filteredProbesCalc = function(covMatrixRow, rd){
  CpGRD1 = sum(covMatrixRow[pheno == levels(pheno)[1]] > rd, na.rm = TRUE)
  CpGRD2 = sum(covMatrixRow[pheno == levels(pheno)[2]] > rd, na.rm = TRUE)
  CpGRD = ifelse(CpGRD1 > nSampleNeeded & CpGRD2 > nSampleNeeded, T, F)
  return(CpGRD)
}

ni = 50
benchmark = function(subSamples){
  nSmol = nTime = rSmol = rTime = matrix(nrow = ni, ncol = 1)
  for (i in 1:ni){
    subsetMatrix = RRBSMatrix[sample(1:nrow(RRBSMatrix),subSamples),covCols]

    x = Sys.time()
    CpGRD = apply(subsetMatrix, 1, filteredProbesCalc, rd)
    nSmol[i] = sum(CpGRD > 0)/nrow(subsetMatrix)
    nTime[i] = Sys.time() - x

    x = Sys.time()
    rSmol[i] = parameterChecksFunction(subsetMatrix)[4]
    rTime[i] = Sys.time() - x
  }
  message(paste("working on", subSamples, "sites"))
  return(c(nSmol, nTime, rSmol, rTime))
}

subSamples = c(1000, seq(2500, 50000, 2500))
benchmarkResults = sapply(subSamples, benchmark)

save(benchmarkResults, file = "/mnt/data1/Thea/Simulations/data/benchmarkingDataAll.Rdata")


